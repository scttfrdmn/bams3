package bams3

import (
	"crypto/sha256"
	"encoding/json"
	"fmt"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"
	"sync"
	"time"
)

// ParallelWriter writes BAMS3 datasets with parallel chunk compression
type ParallelWriter struct {
	*Writer                       // Embed base writer
	workers        int            // Number of parallel workers
	chunkQueue     chan chunkJob  // Queue of chunks to compress
	resultQueue    chan chunkResult // Queue of compressed results
	writerWg       sync.WaitGroup // Wait group for workers
	resultWg       sync.WaitGroup // Wait group for result writer
	errorChan      chan error     // Error channel
	stopOnce       sync.Once      // Ensure cleanup happens once
	chunkResults   []chunkResult  // Accumulated results for metadata
	resultsMutex   sync.Mutex     // Protect chunkResults
}

// chunkJob represents a chunk compression job
type chunkJob struct {
	key   string      // Chunk key (e.g., "chr1:0:1048576")
	reads []ChunkRead // Reads in this chunk
	index int         // Original index for ordering
}

// chunkResult represents a compressed chunk result
type chunkResult struct {
	key   string    // Chunk key
	reads []ChunkRead // Original reads (for writing)
	data  []byte    // Compressed chunk data
	index int       // Original index for ordering
	err   error     // Error if compression failed
}

// NewParallelWriter creates a new parallel BAMS3 writer
func NewParallelWriter(path string, chunkSize int, compression string, format string, workers int) (*ParallelWriter, error) {
	// Default to number of CPUs if workers not specified
	if workers <= 0 {
		workers = runtime.NumCPU()
	}

	// Limit workers to avoid memory exhaustion
	maxWorkers := 32
	if workers > maxWorkers {
		workers = maxWorkers
	}

	// Create base writer
	baseWriter, err := NewWriter(path, chunkSize, compression, format)
	if err != nil {
		return nil, err
	}

	pw := &ParallelWriter{
		Writer:       baseWriter,
		workers:      workers,
		chunkQueue:   make(chan chunkJob, workers*2), // Buffer for smooth pipeline
		resultQueue:  make(chan chunkResult, workers*2),
		errorChan:    make(chan error, 1), // Buffered for non-blocking send
		chunkResults: make([]chunkResult, 0, 100),
	}

	return pw, nil
}

// Start starts the worker pool
func (pw *ParallelWriter) Start() {
	// Start compression workers
	for i := 0; i < pw.workers; i++ {
		pw.writerWg.Add(1)
		go pw.compressionWorker(i)
	}

	// Start result writer (collects compressed chunks)
	pw.resultWg.Add(1)
	go pw.resultCollector()
}

// compressionWorker compresses chunks in parallel
func (pw *ParallelWriter) compressionWorker(id int) {
	defer pw.writerWg.Done()

	for job := range pw.chunkQueue {
		// Check for early termination
		select {
		case err := <-pw.errorChan:
			// Error already reported, just exit
			_ = err
			return
		default:
		}

		// Compress chunk
		result := chunkResult{
			key:   job.key,
			reads: job.reads,
			index: job.index,
		}

		// Encode chunk based on format
		var data []byte
		var err error

		if pw.format == "binary" && pw.binaryWriter != nil {
			// Binary format (compression handled by binary writer)
			data, err = pw.binaryWriter.WriteChunk(job.reads)
		} else {
			// JSON format
			data, err = json.Marshal(job.reads)
			if err == nil && pw.compressor != nil && pw.compression == "zstd" {
				// Compress JSON
				data, err = pw.compressor.Compress(data)
			}
		}

		if err != nil {
			result.err = fmt.Errorf("worker %d failed to compress chunk %s: %w", id, job.key, err)
			// Send error result
			select {
			case pw.resultQueue <- result:
			case <-pw.errorChan:
				return
			}
			continue
		}

		result.data = data

		// Send successful result
		select {
		case pw.resultQueue <- result:
		case <-pw.errorChan:
			return
		}
	}
}

// resultCollector collects compressed chunks (order preservation for metadata)
func (pw *ParallelWriter) resultCollector() {
	defer pw.resultWg.Done()

	for result := range pw.resultQueue {
		if result.err != nil {
			// Report error and signal stop
			select {
			case pw.errorChan <- result.err:
			default:
			}
			return
		}

		// Accumulate results
		pw.resultsMutex.Lock()
		pw.chunkResults = append(pw.chunkResults, result)
		pw.resultsMutex.Unlock()
	}
}

// SubmitChunk submits a chunk for parallel compression
func (pw *ParallelWriter) SubmitChunk(key string, reads []ChunkRead, index int) {
	pw.chunkQueue <- chunkJob{
		key:   key,
		reads: reads,
		index: index,
	}
}

// FinalizeParallel waits for all workers to complete and writes chunks to disk
func (pw *ParallelWriter) FinalizeParallel() error {
	// Close chunk queue (no more jobs)
	close(pw.chunkQueue)

	// Wait for all compression workers to finish
	pw.writerWg.Wait()

	// Check for errors from workers
	select {
	case err := <-pw.errorChan:
		close(pw.resultQueue)
		return err
	default:
	}

	// Close result queue
	close(pw.resultQueue)

	// Wait for result collector to finish
	pw.resultWg.Wait()

	// Check for errors from collector
	select {
	case err := <-pw.errorChan:
		return err
	default:
	}

	// Now write all chunks to disk sequentially (to generate correct metadata)
	fmt.Printf("Writing %d chunks to disk...\n", len(pw.chunkResults))

	for _, result := range pw.chunkResults {
		if err := pw.writeChunkToDisk(result.key, result.reads, result.data); err != nil {
			return fmt.Errorf("failed to write chunk %s: %w", result.key, err)
		}
	}

	// Write header
	if err := pw.writeHeader(); err != nil {
		return fmt.Errorf("failed to write header: %w", err)
	}

	// Calculate mean coverage
	totalRefLength := int64(0)
	for _, sq := range pw.header.SQ {
		if lnStr, ok := sq["LN"]; ok {
			var ln int64
			fmt.Sscanf(lnStr, "%d", &ln)
			totalRefLength += ln
		}
	}
	if totalRefLength > 0 {
		pw.statistics.MeanCoverage = float64(pw.statistics.TotalBases) / float64(totalRefLength)
	}

	// Write metadata
	pw.metadata.Statistics = pw.statistics
	pw.metadata.Chunks = pw.chunkInfos
	if err := pw.writeMetadata(); err != nil {
		return fmt.Errorf("failed to write metadata: %w", err)
	}

	// Write spatial index
	if err := pw.writeSpatialIndex(); err != nil {
		return fmt.Errorf("failed to write spatial index: %w", err)
	}

	// Close compressor if used
	if pw.compressor != nil {
		pw.compressor.Close()
	}

	fmt.Println("✓ Conversion complete!")
	fmt.Printf("\nDataset summary:\n")
	fmt.Printf("  Location: %s\n", pw.path)
	fmt.Printf("  Total reads: %d\n", pw.statistics.TotalReads)
	fmt.Printf("  Mapped reads: %d\n", pw.statistics.MappedReads)
	fmt.Printf("  Chunks: %d\n", len(pw.chunkInfos))
	fmt.Printf("  Workers: %d\n", pw.workers)
	if pw.compression != "none" {
		fmt.Printf("  Compression: %s\n", pw.compression)
	}

	return nil
}

// writeChunkToDisk writes a pre-compressed chunk to disk and generates metadata
func (pw *ParallelWriter) writeChunkToDisk(chunkKey string, reads []ChunkRead, compressedData []byte) error {
	// Parse chunk key
	var reference string
	var start, end int

	if chunkKey == "unmapped:0:0" {
		reference = "unmapped"
		start = 0
		end = 0
	} else {
		parts := strings.Split(chunkKey, ":")
		if len(parts) != 3 {
			return fmt.Errorf("invalid chunk key format: %s", chunkKey)
		}
		reference = parts[0]
		var err error
		start, err = strconv.Atoi(parts[1])
		if err != nil {
			return fmt.Errorf("failed to parse start position: %w", err)
		}
		end, err = strconv.Atoi(parts[2])
		if err != nil {
			return fmt.Errorf("failed to parse end position: %w", err)
		}
	}

	// Determine chunk path
	var chunkPath string
	if reference == "unmapped" {
		chunkPath = filepath.Join("data", "unmapped.chunk")
	} else {
		refDir := filepath.Join("data", reference)
		if err := pw.storage.MkdirAll(refDir); err != nil {
			return err
		}
		chunkPath = filepath.Join(refDir, fmt.Sprintf("%09d-%09d.chunk", start, end))
	}

	// Write compressed data to file
	if err := pw.storage.WriteFile(chunkPath, compressedData); err != nil {
		return err
	}

	// Calculate checksum (of compressed data)
	hash := sha256.Sum256(compressedData)
	checksum := fmt.Sprintf("%x", hash)

	// Get file size
	fileSize := int64(len(compressedData))

	// Record chunk info
	chunkInfo := ChunkInfo{
		Path:        filepath.Join("data", reference, filepath.Base(chunkPath)),
		Reference:   reference,
		Start:       start,
		End:         end,
		Reads:       len(reads),
		SizeBytes:   fileSize,
		Compression: pw.compression,
		Checksum:    checksum,
		Created:     time.Now(),
	}

	if reference == "unmapped" {
		chunkInfo.Path = filepath.Join("data", "unmapped.chunk")
	}

	pw.chunkInfos = append(pw.chunkInfos, chunkInfo)

	return nil
}

// GetWorkerCount returns the number of parallel workers
func (pw *ParallelWriter) GetWorkerCount() int {
	return pw.workers
}
