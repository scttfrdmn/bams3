package bams3

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"sync"
	"time"
)

// StreamConfig holds configuration for streaming conversion
type StreamConfig struct {
	// Resource allocation
	Workers         int   // Number of parallel workers (default: runtime.NumCPU())
	BufferSize      int64 // Read buffer size in bytes (default: 25% of RAM)
	SortBufferSize  int64 // In-memory sort buffer in bytes (default: min(8GB, 25% RAM))
	ChunkSize       int   // Chunk size in base pairs (default: 1MB)

	// Compression
	Compression      string // Compression algorithm (default: "zstd")
	CompressionLevel int    // Compression level (default: 3)
	Format          string  // Output format (default: "binary")

	// Disk spill (for large datasets)
	SpillToDisk bool   // Enable disk spill for sorting (default: true)
	TempDir     string // Temporary directory for spills (default: os.TempDir())

	// S3 optimization
	S3PartSize    int64 // S3 multipart upload part size (default: 10MB)
	S3Concurrency int   // S3 upload concurrency (default: Workers * 2)

	// Progress reporting
	ShowProgress bool          // Show progress bar (default: true)
	ProgressInterval time.Duration // Progress update interval (default: 1s)

	// Advanced options
	Prefetch bool // Prefetch next chunk metadata (default: true)

	// Computed fields (not user-configurable)
	availableMemory int64
}

// NewStreamConfig creates a StreamConfig with smart defaults
func NewStreamConfig() *StreamConfig {
	memStats := getSystemMemory()

	workers := detectOptimalWorkers() // Use P-cores on hybrid architectures
	bufferSize := memStats.Total / 4  // 25% of RAM
	sortBufferSize := min64(8*GB, memStats.Total/4)

	return &StreamConfig{
		Workers:          workers,
		BufferSize:       bufferSize,
		SortBufferSize:   sortBufferSize,
		ChunkSize:        1048576, // 1MB
		Compression:      "zstd",
		CompressionLevel: 3,
		Format:           "binary",
		SpillToDisk:      true,
		TempDir:          os.TempDir(),
		S3PartSize:       10 * MB,
		S3Concurrency:    workers * 2,
		ShowProgress:     true,
		ProgressInterval: 1 * time.Second,
		Prefetch:         true,
		availableMemory:  memStats.Available,
	}
}

// Validate checks configuration and warns about potential issues
func (c *StreamConfig) Validate() error {
	if c.Workers < 1 {
		return fmt.Errorf("workers must be >= 1")
	}
	if c.Workers > 64 {
		fmt.Fprintf(os.Stderr, "Warning: Workers > 64 may cause diminishing returns\n")
	}
	if c.SortBufferSize < 100*MB {
		fmt.Fprintf(os.Stderr, "Warning: Sort buffer < 100MB may cause excessive disk spills\n")
	}
	if c.BufferSize > c.availableMemory {
		return fmt.Errorf("buffer size (%d GB) exceeds available memory (%d GB)",
			c.BufferSize/GB, c.availableMemory/GB)
	}
	if c.ChunkSize < 64*1024 || c.ChunkSize > 16*1024*1024 {
		return fmt.Errorf("chunk size must be between 64KB and 16MB")
	}
	return nil
}

// ShowConfig prints the effective configuration
func (c *StreamConfig) ShowConfig() {
	// Get system memory info for context
	memStats := getSystemMemory()

	fmt.Fprintf(os.Stderr, "System Information:\n")
	fmt.Fprintf(os.Stderr, "  Total RAM: %.1f GB\n", float64(memStats.Total)/float64(GB))
	fmt.Fprintf(os.Stderr, "  Available RAM: %.1f GB\n", float64(memStats.Available)/float64(GB))

	// Show CPU core information
	totalCores := runtime.NumCPU()
	optimalWorkers := detectOptimalWorkers()
	if optimalWorkers < totalCores {
		fmt.Fprintf(os.Stderr, "  CPU cores: %d total (%d performance, %d efficiency)\n",
			totalCores, optimalWorkers, totalCores-optimalWorkers)
	} else {
		fmt.Fprintf(os.Stderr, "  CPU cores: %d\n", totalCores)
	}
	fmt.Fprintf(os.Stderr, "\n")

	fmt.Fprintf(os.Stderr, "Configuration:\n")
	fmt.Fprintf(os.Stderr, "  Workers: %d\n", c.Workers)
	fmt.Fprintf(os.Stderr, "  Buffer: %.1f GB\n", float64(c.BufferSize)/float64(GB))
	fmt.Fprintf(os.Stderr, "  Sort buffer: %.1f GB\n", float64(c.SortBufferSize)/float64(GB))
	fmt.Fprintf(os.Stderr, "  Chunk size: %s\n", FormatChunkSize(c.ChunkSize))
	fmt.Fprintf(os.Stderr, "  Compression: %s (level %d)\n", c.Compression, c.CompressionLevel)
	fmt.Fprintf(os.Stderr, "  Format: %s\n", c.Format)
	if c.SpillToDisk {
		fmt.Fprintf(os.Stderr, "  Disk spill: enabled (%s)\n", c.TempDir)
	} else {
		fmt.Fprintf(os.Stderr, "  Disk spill: disabled\n")
	}
	fmt.Fprintf(os.Stderr, "\n")
}

// SpillFile tracks a temporary sorted file on disk
type SpillFile struct {
	Path      string
	NumReads  int
	SizeBytes int64
}

// StreamConverter handles streaming conversion from stdin to BAMS3
type StreamConverter struct {
	config     *StreamConfig
	writer     *ParallelWriter
	stats      ConversionStats
	statsMutex sync.Mutex

	// Read frontier tracking for incremental flushing
	// Maps reference ID to the maximum position processed
	readFrontier map[int]int

	// Disk spill tracking for external sort
	spillFiles []SpillFile
	spillDir   string
}

// ConversionStats tracks conversion progress
type ConversionStats struct {
	ReadsProcessed  int64
	ReadsInSort     int64
	ChunksFinalized int64
	ChunksUploaded  int64
	BytesRead       int64
	BytesWritten    int64
	SpillsToDisK    int64
	StartTime       time.Time
	LastUpdate      time.Time
}

// NewStreamConverter creates a new streaming converter
func NewStreamConverter(outputPath string, config *StreamConfig) (*StreamConverter, error) {
	if config == nil {
		config = NewStreamConfig()
	}

	if err := config.Validate(); err != nil {
		return nil, err
	}

	// Create parallel writer
	writer, err := NewParallelWriter(
		outputPath,
		config.ChunkSize,
		config.Compression,
		config.Format,
		config.Workers,
	)
	if err != nil {
		return nil, fmt.Errorf("failed to create writer: %w", err)
	}

	// Create temporary directory for disk spills if enabled
	var spillDir string
	if config.SpillToDisk {
		spillDir, err = os.MkdirTemp(config.TempDir, "bams3-spill-*")
		if err != nil {
			return nil, fmt.Errorf("failed to create spill directory: %w", err)
		}
	}

	return &StreamConverter{
		config:       config,
		writer:       writer,
		readFrontier: make(map[int]int),
		spillFiles:   make([]SpillFile, 0),
		spillDir:     spillDir,
		stats: ConversionStats{
			StartTime:  time.Now(),
			LastUpdate: time.Now(),
		},
	}, nil
}

// ConvertFromStream converts SAM/BAM from stdin to BAMS3
func (sc *StreamConverter) ConvertFromStream(reader io.Reader, source Source) error {
	fmt.Fprintf(os.Stderr, "Converting from stdin...\n")
	sc.config.ShowConfig()

	// Ensure cleanup of temporary files
	defer sc.cleanup()

	// Parse header from input (returns scanner positioned at first data line)
	header, refNameToID, scanner, err := sc.parseSAMHeader(reader)
	if err != nil {
		return fmt.Errorf("failed to parse header: %w", err)
	}

	// Set header and source
	sc.writer.SetHeader(header)
	sc.writer.SetSource(source)

	// Start parallel writer workers
	sc.writer.Start()

	// Start progress reporter
	progressDone := make(chan bool)
	if sc.config.ShowProgress {
		go sc.reportProgress(progressDone)
	}

	// Read and process records (continue from where header parsing left off)
	if err := sc.readAndSort(scanner, refNameToID); err != nil {
		close(progressDone)
		return fmt.Errorf("failed to read and sort: %w", err)
	}

	// Stop progress reporter
	if sc.config.ShowProgress {
		close(progressDone)
	}

	fmt.Fprintf(os.Stderr, "\nTotal reads processed: %d\n", sc.stats.ReadsProcessed)
	fmt.Fprintf(os.Stderr, "Finalizing chunks and uploading...\n\n")

	// Finalize parallel writer
	if err := sc.writer.FinalizeParallel(); err != nil {
		return fmt.Errorf("failed to finalize: %w", err)
	}

	// Print final stats
	sc.printFinalStats()

	return nil
}

// parseSAMHeader parses SAM header lines from input
// Returns header, refName mapping, scanner positioned at first data line, and error
func (sc *StreamConverter) parseSAMHeader(reader io.Reader) (Header, map[string]int, *bufio.Scanner, error) {
	header := Header{
		HD: make(map[string]string),
		SQ: []map[string]string{},
		RG: []map[string]string{},
		PG: []map[string]string{},
		CO: []string{},
	}

	refNameToID := make(map[string]int)
	refID := 0

	scanner := bufio.NewScanner(reader)
	scanner.Buffer(make([]byte, 0, int(sc.config.BufferSize)), int(sc.config.BufferSize))

	for scanner.Scan() {
		line := scanner.Text()

		// Stop at first non-header line (but scanner is now positioned at this line)
		if !strings.HasPrefix(line, "@") {
			// Process this line in readAndSort
			break
		}

		fields := strings.Split(line, "\t")
		if len(fields) < 1 {
			continue
		}

		tag := fields[0]
		switch tag {
		case "@HD":
			for _, field := range fields[1:] {
				parts := strings.SplitN(field, ":", 2)
				if len(parts) == 2 {
					header.HD[parts[0]] = parts[1]
				}
			}

		case "@SQ":
			sq := make(map[string]string)
			var refName string
			for _, field := range fields[1:] {
				parts := strings.SplitN(field, ":", 2)
				if len(parts) == 2 {
					sq[parts[0]] = parts[1]
					if parts[0] == "SN" {
						refName = parts[1]
					}
				}
			}
			header.SQ = append(header.SQ, sq)
			if refName != "" {
				refNameToID[refName] = refID
				refID++
			}

		case "@RG":
			rg := make(map[string]string)
			for _, field := range fields[1:] {
				parts := strings.SplitN(field, ":", 2)
				if len(parts) == 2 {
					rg[parts[0]] = parts[1]
				}
			}
			header.RG = append(header.RG, rg)

		case "@PG":
			pg := make(map[string]string)
			for _, field := range fields[1:] {
				parts := strings.SplitN(field, ":", 2)
				if len(parts) == 2 {
					pg[parts[0]] = parts[1]
				}
			}
			header.PG = append(header.PG, pg)

		case "@CO":
			if len(fields) > 1 {
				header.CO = append(header.CO, strings.Join(fields[1:], "\t"))
			}
		}
	}

	if err := scanner.Err(); err != nil {
		return header, nil, nil, err
	}

	// Return scanner so readAndSort can continue from current position
	return header, refNameToID, scanner, nil
}

// readAndSort reads records from scanner and sorts them
func (sc *StreamConverter) readAndSort(scanner *bufio.Scanner, refNameToID map[string]int) error {
	var sortBuffer []Read
	sortBufferBytes := int64(0)

	// Process current line first (already scanned by header parser)
	line := scanner.Text()
	if !strings.HasPrefix(line, "@") {
		read, err := parseSAMLine(line, refNameToID)
		if err == nil {
			sc.updateStats(1, int64(len(line)))
			sortBuffer = append(sortBuffer, read)
			sortBufferBytes += int64(len(read.Sequence) + len(read.Quality) + 100)
		}
	}

	// Continue with remaining lines
	for scanner.Scan() {
		line := scanner.Text()

		// Skip any remaining header lines (safety check)
		if strings.HasPrefix(line, "@") {
			continue
		}

		// Parse SAM record
		read, err := parseSAMLine(line, refNameToID)
		if err != nil {
			// Skip malformed lines with warning
			fmt.Fprintf(os.Stderr, "Warning: skipping malformed line: %v\n", err)
			continue
		}

		sc.updateStats(1, int64(len(line)))

		// Add to sort buffer
		sortBuffer = append(sortBuffer, read)
		sortBufferBytes += int64(len(read.Sequence) + len(read.Quality) + 100) // Rough estimate

		// Check if sort buffer is full
		if sortBufferBytes >= sc.config.SortBufferSize {
			if err := sc.flushSortBuffer(&sortBuffer, &sortBufferBytes, false); err != nil {
				return err
			}
		}
	}

	if err := scanner.Err(); err != nil {
		return fmt.Errorf("error reading input: %w", err)
	}

	// Final merge and flush
	if err := sc.finalMergeAndFlush(&sortBuffer); err != nil {
		return err
	}

	return nil
}

// flushSortBuffer sorts and writes records from the sort buffer
// If buffer is full but chunks can't be flushed (due to frontier), spills to disk
func (sc *StreamConverter) flushSortBuffer(buffer *[]Read, bufferBytes *int64, isFinal bool) error {
	if len(*buffer) == 0 {
		return nil
	}

	// Sort by coordinate
	sort.Slice(*buffer, func(i, j int) bool {
		if (*buffer)[i].ReferenceID != (*buffer)[j].ReferenceID {
			return (*buffer)[i].ReferenceID < (*buffer)[j].ReferenceID
		}
		return (*buffer)[i].Position < (*buffer)[j].Position
	})

	sc.statsMutex.Lock()
	sc.stats.ReadsInSort = int64(len(*buffer))
	sc.statsMutex.Unlock()

	// Update read frontier (track maximum position per reference in this batch)
	// This tells us how far we've progressed on each reference
	batchFrontier := make(map[int]int)
	for i := range *buffer {
		refID := (*buffer)[i].ReferenceID
		pos := (*buffer)[i].Position
		if currentMax, exists := batchFrontier[refID]; !exists || pos > currentMax {
			batchFrontier[refID] = pos
		}
	}

	// Update global frontier (maximum position we've seen across all batches)
	for refID, maxPos := range batchFrontier {
		if currentFrontier, exists := sc.readFrontier[refID]; !exists || maxPos > currentFrontier {
			sc.readFrontier[refID] = maxPos
		}
	}

	initialBufferSize := len(*buffer)

	// First, check if we can do incremental flush by examining existing chunks
	safetyMargin := sc.config.ChunkSize * 2
	existingChunks := sc.writer.Writer.GetChunks()
	canFlushExisting := false

	for chunkKey := range existingChunks {
		if sc.shouldFlushChunk(chunkKey, safetyMargin) {
			canFlushExisting = true
			break
		}
	}

	// If we can't flush existing chunks and buffer is full, we need to spill
	// (unless this is preventing us from making progress)
	if !canFlushExisting && len(existingChunks) > 0 && sc.config.SpillToDisk && sc.spillDir != "" {
		// Spill current buffer to disk without adding to writer
		fmt.Fprintf(os.Stderr, "\nSpilling %d reads to disk (spill #%d)...\n",
			initialBufferSize, len(sc.spillFiles)+1)

		spillNum := len(sc.spillFiles)
		spillWriter, err := NewSpillWriter(sc.spillDir, spillNum)
		if err != nil {
			return fmt.Errorf("failed to create spill writer: %w", err)
		}

		for i := range *buffer {
			if err := spillWriter.WriteRead(&(*buffer)[i]); err != nil {
				spillWriter.Close()
				return fmt.Errorf("failed to write read to spill: %w", err)
			}
		}

		spillFile, err := spillWriter.Close()
		if err != nil {
			return fmt.Errorf("failed to close spill writer: %w", err)
		}

		sc.spillFiles = append(sc.spillFiles, spillFile)

		sc.statsMutex.Lock()
		sc.stats.SpillsToDisK++
		sc.statsMutex.Unlock()

		fmt.Fprintf(os.Stderr, "Spilled %d reads (%.1f MB) to %s\n",
			spillFile.NumReads, float64(spillFile.SizeBytes)/float64(MB), spillFile.Path)

		// Clear buffer - reads are now on disk
		*buffer = (*buffer)[:0]
		*bufferBytes = 0

		return nil
	}

	// Not spilling - add reads to writer for incremental flush
	for i := range *buffer {
		refName := sc.getRefName((*buffer)[i].ReferenceID)
		if err := sc.writer.AddRead((*buffer)[i], refName); err != nil {
			return fmt.Errorf("failed to add read: %w", err)
		}
	}

	// Try to flush chunks that are behind the frontier
	chunks := sc.writer.Writer.GetChunks()
	chunkIndex := 0
	flushedAny := false

	for chunkKey, reads := range chunks {
		if sc.shouldFlushChunk(chunkKey, safetyMargin) {
			sc.writer.SubmitChunk(chunkKey, reads, chunkIndex)
			chunkIndex++
			flushedAny = true

			sc.statsMutex.Lock()
			sc.stats.ChunksFinalized++
			sc.statsMutex.Unlock()
		}
	}

	// Clear buffer regardless (reads are now in writer's chunks)
	*buffer = (*buffer)[:0]
	*bufferBytes = 0

	if !flushedAny && len(existingChunks) > 10 {
		// Warning: accumulating many chunks without flushing
		fmt.Fprintf(os.Stderr, "\nWarning: %d chunks accumulated without flushing\n", len(chunks))
	}

	return nil
}

// finalMergeAndFlush performs final merge of spill files and in-memory buffer
func (sc *StreamConverter) finalMergeAndFlush(buffer *[]Read) error {
	// If no spill files, just do a normal final flush
	if len(sc.spillFiles) == 0 {
		bufferBytes := int64(0)
		for i := range *buffer {
			bufferBytes += int64(len((*buffer)[i].Sequence) + len((*buffer)[i].Quality) + 100)
		}
		return sc.flushSortBuffer(buffer, &bufferBytes, true)
	}

	// We have spill files - need to do k-way merge
	fmt.Fprintf(os.Stderr, "\nPerforming k-way merge of %d spill files + in-memory buffer...\n",
		len(sc.spillFiles))

	// Open all spill readers
	spillReaders := make([]*SpillReader, 0, len(sc.spillFiles))
	for i, spillFile := range sc.spillFiles {
		reader, err := NewSpillReader(spillFile.Path)
		if err != nil {
			// Close any readers we've already opened
			for _, r := range spillReaders {
				r.Close()
			}
			return fmt.Errorf("failed to open spill file %d: %w", i, err)
		}
		spillReaders = append(spillReaders, reader)
	}

	// Ensure readers are closed when done
	defer func() {
		for _, reader := range spillReaders {
			reader.Close()
		}
	}()

	// Sort in-memory buffer first
	sort.Slice(*buffer, func(i, j int) bool {
		if (*buffer)[i].ReferenceID != (*buffer)[j].ReferenceID {
			return (*buffer)[i].ReferenceID < (*buffer)[j].ReferenceID
		}
		return (*buffer)[i].Position < (*buffer)[j].Position
	})

	// Perform k-way merge
	mergedReads, err := kWayMerge(spillReaders, *buffer)
	if err != nil {
		return fmt.Errorf("k-way merge failed: %w", err)
	}

	fmt.Fprintf(os.Stderr, "Merged %d total reads, writing to chunks...\n", len(mergedReads))

	// Write all merged reads to chunks
	for i := range mergedReads {
		refName := sc.getRefName(mergedReads[i].ReferenceID)
		if err := sc.writer.AddRead(mergedReads[i], refName); err != nil {
			return fmt.Errorf("failed to add merged read: %w", err)
		}
	}

	// Flush all accumulated chunks
	chunks := sc.writer.Writer.GetChunks()
	chunkIndex := 0

	for chunkKey, reads := range chunks {
		sc.writer.SubmitChunk(chunkKey, reads, chunkIndex)
		chunkIndex++

		sc.statsMutex.Lock()
		sc.stats.ChunksFinalized++
		sc.statsMutex.Unlock()
	}

	// Clear buffer
	*buffer = (*buffer)[:0]

	return nil
}

// shouldFlushChunk determines if a chunk is safe to flush based on frontier
func (sc *StreamConverter) shouldFlushChunk(chunkKey string, safetyMargin int) bool {
	// Parse chunk key format: "refName:start-end"
	parts := strings.Split(chunkKey, ":")
	if len(parts) != 2 {
		return true // Malformed key, flush it
	}

	refName := parts[0]
	positions := strings.Split(parts[1], "-")
	if len(positions) != 2 {
		return true
	}

	// Parse end position
	endPos, err := strconv.Atoi(positions[1])
	if err != nil {
		return true
	}

	// Get reference ID for this chunk
	refID := sc.getRefIDFromName(refName)
	if refID < 0 {
		return true // Unmapped or unknown reference
	}

	// Check if this chunk is safely behind the frontier
	if frontier, exists := sc.readFrontier[refID]; exists {
		// Flush if chunk end is well behind our current position (frontier - safety margin)
		// frontier = max position we've processed
		// If chunk ends before (frontier - margin), it's safe to flush
		safeFlushPoint := frontier - safetyMargin
		if endPos <= safeFlushPoint {
			return true
		}
		// Chunk is too close to or ahead of our current position, keep it
		return false
	}

	// No frontier data for this reference yet
	// This means we haven't processed any reads for this reference
	// Conservative: don't flush until we have frontier info
	return false
}

// getRefIDFromName returns reference ID for reference name
func (sc *StreamConverter) getRefIDFromName(refName string) int {
	if refName == "*" {
		return -1
	}

	// Search through SQ records
	for i, sq := range sc.writer.Writer.header.SQ {
		if sq["SN"] == refName {
			return i
		}
	}

	return -1
}

// getRefName returns reference name for reference ID
func (sc *StreamConverter) getRefName(refID int) string {
	if refID < 0 {
		return "*"
	}

	// Get from header
	if refID < len(sc.writer.Writer.header.SQ) {
		return sc.writer.Writer.header.SQ[refID]["SN"]
	}

	return "*"
}

// updateStats updates conversion statistics
func (sc *StreamConverter) updateStats(reads int64, bytes int64) {
	sc.statsMutex.Lock()
	defer sc.statsMutex.Unlock()

	sc.stats.ReadsProcessed += reads
	sc.stats.BytesRead += bytes
}

// reportProgress prints progress updates
func (sc *StreamConverter) reportProgress(done chan bool) {
	ticker := time.NewTicker(sc.config.ProgressInterval)
	defer ticker.Stop()

	for {
		select {
		case <-done:
			return
		case <-ticker.C:
			sc.printProgress()
		}
	}
}

// printProgress prints current progress
func (sc *StreamConverter) printProgress() {
	sc.statsMutex.Lock()
	defer sc.statsMutex.Unlock()

	elapsed := time.Since(sc.stats.StartTime)
	readsPerSec := float64(sc.stats.ReadsProcessed) / elapsed.Seconds()
	mbPerSec := float64(sc.stats.BytesRead) / elapsed.Seconds() / MB

	fmt.Fprintf(os.Stderr, "\rProgress: %d reads (%.1f K reads/s, %.1f MB/s) | Chunks: %d finalized | Elapsed: %s",
		sc.stats.ReadsProcessed,
		readsPerSec/1000,
		mbPerSec,
		sc.stats.ChunksFinalized,
		formatDuration(elapsed),
	)
}

// printFinalStats prints final conversion statistics
func (sc *StreamConverter) printFinalStats() {
	elapsed := time.Since(sc.stats.StartTime)
	readsPerSec := float64(sc.stats.ReadsProcessed) / elapsed.Seconds()

	fmt.Fprintf(os.Stderr, "\n")
	fmt.Fprintf(os.Stderr, "Conversion complete!\n")
	fmt.Fprintf(os.Stderr, "  Total reads: %d\n", sc.stats.ReadsProcessed)
	fmt.Fprintf(os.Stderr, "  Total chunks: %d\n", sc.stats.ChunksFinalized)
	fmt.Fprintf(os.Stderr, "  Throughput: %.1f K reads/s\n", readsPerSec/1000)
	fmt.Fprintf(os.Stderr, "  Elapsed time: %s\n", formatDuration(elapsed))
	if sc.stats.SpillsToDisK > 0 {
		fmt.Fprintf(os.Stderr, "  Disk spills: %d\n", sc.stats.SpillsToDisK)
	}
}

// parseSAMLine parses a SAM format line into a Read
func parseSAMLine(line string, refNameToID map[string]int) (Read, error) {
	fields := strings.Split(line, "\t")
	if len(fields) < 11 {
		return Read{}, fmt.Errorf("invalid SAM line: too few fields")
	}

	read := Read{
		Name:     fields[0],
		Sequence: fields[9],
		Quality:  fields[10],
		Tags:     make(map[string]interface{}),
	}

	// Parse flag
	fmt.Sscanf(fields[1], "%d", &read.Flag)

	// Parse reference ID from RNAME (fields[2])
	rname := fields[2]
	if rname == "*" {
		read.ReferenceID = -1
	} else {
		if refID, ok := refNameToID[rname]; ok {
			read.ReferenceID = refID
		} else {
			// Unknown reference - treat as unmapped
			read.ReferenceID = -1
		}
	}

	// Parse position
	fmt.Sscanf(fields[3], "%d", &read.Position)

	// Parse mapping quality
	var mapq int
	fmt.Sscanf(fields[4], "%d", &mapq)
	read.MappingQuality = uint8(mapq)

	// Parse CIGAR
	read.CIGAR = fields[5]

	// Parse mate reference (RNEXT, fields[6])
	rnext := fields[6]
	if rnext == "*" {
		read.MateReferenceID = -1
	} else if rnext == "=" {
		read.MateReferenceID = read.ReferenceID
	} else {
		if refID, ok := refNameToID[rnext]; ok {
			read.MateReferenceID = refID
		} else {
			read.MateReferenceID = -1
		}
	}

	// Parse mate position (PNEXT, fields[7])
	fmt.Sscanf(fields[7], "%d", &read.MatePosition)

	// Parse template length (TLEN, fields[8])
	fmt.Sscanf(fields[8], "%d", &read.TemplateLength)

	// Parse tags (fields 11+)
	for i := 11; i < len(fields); i++ {
		parts := strings.SplitN(fields[i], ":", 3)
		if len(parts) == 3 {
			read.Tags[parts[0]] = parts[2]
		}
	}

	return read, nil
}

// Helper functions

func min64(a, b int64) int64 {
	if a < b {
		return a
	}
	return b
}

func formatDuration(d time.Duration) string {
	d = d.Round(time.Second)
	h := d / time.Hour
	d -= h * time.Hour
	m := d / time.Minute
	d -= m * time.Minute
	s := d / time.Second

	if h > 0 {
		return fmt.Sprintf("%dh %dm %ds", h, m, s)
	}
	if m > 0 {
		return fmt.Sprintf("%dm %ds", m, s)
	}
	return fmt.Sprintf("%ds", s)
}

// Constants
const (
	KB = 1024
	MB = 1024 * KB
	GB = 1024 * MB
)

// SystemMemory holds system memory information
type SystemMemory struct {
	Total     int64
	Available int64
	Used      int64
}

// getSystemMemory returns system memory stats
func getSystemMemory() SystemMemory {
	total, available := detectSystemMemory()

	// Fallback to sensible defaults if detection fails
	if total == 0 {
		total = 16 * GB
		available = 12 * GB
	}

	used := total - available
	return SystemMemory{
		Total:     total,
		Available: available,
		Used:      used,
	}
}

// cleanup removes temporary spill files and directory
func (sc *StreamConverter) cleanup() {
	if sc.spillDir == "" {
		return
	}

	// Remove all spill files
	for _, spill := range sc.spillFiles {
		os.Remove(spill.Path)
	}

	// Remove spill directory
	os.RemoveAll(sc.spillDir)
}
