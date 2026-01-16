package bams3

import (
	"crypto/sha256"
	"encoding/json"
	"fmt"
	"path/filepath"
	"strconv"
	"strings"
	"time"
)

// Writer writes BAMS3 datasets
type Writer struct {
	path          string
	storage       Storage // Storage backend (local or S3)
	chunkSize     int
	compression   string
	format        string // "json" or "binary"
	metadata      Metadata
	header        Header
	chunks        map[string][]ChunkRead // key: "ref:start:end"
	chunkInfos    []ChunkInfo
	statistics    Statistics
	compressor    *Compressor
	binaryWriter  *BinaryChunkWriter
}

// NewWriter creates a new BAMS3 writer
func NewWriter(path string, chunkSize int, compression string, format string) (*Writer, error) {
	// Default format to binary for v0.2
	if format == "" {
		format = "binary"
	}

	// Determine version based on format
	version := "0.1.0"
	if format == "binary" {
		version = "0.2.0"
	}

	// Create storage backend (auto-detects local vs S3)
	storage, err := NewStorage(path)
	if err != nil {
		return nil, fmt.Errorf("failed to create storage backend: %w", err)
	}

	w := &Writer{
		path:        path,
		storage:     storage,
		chunkSize:   chunkSize,
		compression: compression,
		format:      format,
		chunks:      make(map[string][]ChunkRead),
		metadata: Metadata{
			Format:     "bams3",
			Version:    version,
			Created:    time.Now(),
			CreatedBy:  "bams3-go",
			ChunkSize:  chunkSize,
			Compression: CompressionConfig{
				Algorithm: compression,
			},
		},
	}

	// Initialize compressor if compression is enabled
	if compression == "zstd" {
		compressor, err := NewCompressor()
		if err != nil {
			return nil, fmt.Errorf("failed to create compressor: %w", err)
		}
		w.compressor = compressor
	}

	// Initialize binary writer if using binary format
	if format == "binary" {
		binaryWriter, err := NewBinaryChunkWriter(compression)
		if err != nil {
			return nil, fmt.Errorf("failed to create binary writer: %w", err)
		}
		w.binaryWriter = binaryWriter
	}

	// Create directory structure
	if err := storage.MkdirAll(""); err != nil {
		return nil, fmt.Errorf("failed to create dataset directory: %w", err)
	}

	if err := storage.MkdirAll("_index"); err != nil {
		return nil, fmt.Errorf("failed to create index directory: %w", err)
	}

	if err := storage.MkdirAll("data"); err != nil {
		return nil, fmt.Errorf("failed to create data directory: %w", err)
	}

	return w, nil
}

// SetHeader sets the SAM header
func (w *Writer) SetHeader(header Header) {
	w.header = header
}

// SetSource sets the source information
func (w *Writer) SetSource(source Source) {
	w.metadata.Source = source
}

// GetChunks returns the accumulated chunks (for parallel processing)
func (w *Writer) GetChunks() map[string][]ChunkRead {
	return w.chunks
}

// AddRead adds a read to the appropriate chunk
func (w *Writer) AddRead(read Read, refName string) error {
	// Update statistics
	w.statistics.TotalReads++
	if read.ReferenceID >= 0 {
		w.statistics.MappedReads++
	} else {
		w.statistics.UnmappedReads++
	}
	if read.Flag&0x400 != 0 { // Duplicate flag
		w.statistics.DuplicateReads++
	}
	w.statistics.TotalBases += int64(len(read.Sequence))

	// Determine chunk
	var chunkKey string
	if read.ReferenceID < 0 {
		// Unmapped reads
		chunkKey = "unmapped:0:0"
	} else {
		chunkStart := (read.Position / w.chunkSize) * w.chunkSize
		chunkEnd := chunkStart + w.chunkSize
		chunkKey = fmt.Sprintf("%s:%d:%d", refName, chunkStart, chunkEnd)
	}

	// Convert Read to ChunkRead
	chunkRead := ChunkRead{
		Name:            read.Name,
		Flag:            read.Flag,
		ReferenceID:     read.ReferenceID,
		Position:        read.Position,
		MappingQuality:  read.MappingQuality,
		CIGAR:           read.CIGAR,
		Sequence:        read.Sequence,
		Quality:         read.Quality,
		Tags:            read.Tags,
		CIGAROps:        read.CIGAROps, // For binary format
	}

	// Legacy fields for JSON format compatibility
	chunkRead.Ref = read.ReferenceID
	chunkRead.Pos = read.Position
	chunkRead.MapQ = read.MappingQuality
	chunkRead.Seq = read.Sequence
	chunkRead.Qual = read.Quality

	// Add to chunk
	w.chunks[chunkKey] = append(w.chunks[chunkKey], chunkRead)

	return nil
}

// Finalize writes all chunks and metadata
func (w *Writer) Finalize() error {
	fmt.Printf("Writing %d chunks...\n", len(w.chunks))

	// Write each chunk
	for chunkKey, reads := range w.chunks {
		if err := w.writeChunk(chunkKey, reads); err != nil {
			return fmt.Errorf("failed to write chunk %s: %w", chunkKey, err)
		}
	}

	// Write header
	if err := w.writeHeader(); err != nil {
		return fmt.Errorf("failed to write header: %w", err)
	}

	// Calculate mean coverage
	totalRefLength := int64(0)
	for _, sq := range w.header.SQ {
		if lnStr, ok := sq["LN"]; ok {
			var ln int64
			fmt.Sscanf(lnStr, "%d", &ln)
			totalRefLength += ln
		}
	}
	if totalRefLength > 0 {
		w.statistics.MeanCoverage = float64(w.statistics.TotalBases) / float64(totalRefLength)
	}

	// Write metadata
	w.metadata.Statistics = w.statistics
	w.metadata.Chunks = w.chunkInfos
	if err := w.writeMetadata(); err != nil {
		return fmt.Errorf("failed to write metadata: %w", err)
	}

	// Write spatial index
	if err := w.writeSpatialIndex(); err != nil {
		return fmt.Errorf("failed to write spatial index: %w", err)
	}

	// Close compressor if used
	if w.compressor != nil {
		w.compressor.Close()
	}

	fmt.Println("✓ Conversion complete!")
	fmt.Printf("\nDataset summary:\n")
	fmt.Printf("  Location: %s\n", w.path)
	fmt.Printf("  Total reads: %d\n", w.statistics.TotalReads)
	fmt.Printf("  Mapped reads: %d\n", w.statistics.MappedReads)
	fmt.Printf("  Chunks: %d\n", len(w.chunkInfos))
	if w.compression != "none" {
		fmt.Printf("  Compression: %s\n", w.compression)
	}

	return nil
}

// writeChunk writes a chunk to disk
func (w *Writer) writeChunk(chunkKey string, reads []ChunkRead) error {
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
		// Create reference directory
		refDir := filepath.Join("data", reference)
		if err := w.storage.MkdirAll(refDir); err != nil {
			return err
		}
		chunkPath = filepath.Join(refDir, fmt.Sprintf("%09d-%09d.chunk", start, end))
	}

	// Write chunk data - use binary or JSON format
	var data []byte
	compressionType := w.compression

	if w.format == "binary" && w.binaryWriter != nil {
		// Use binary format
		binaryData, err := w.binaryWriter.WriteChunk(reads)
		if err != nil {
			return fmt.Errorf("failed to write binary chunk: %w", err)
		}
		data = binaryData
	} else {
		// Use JSON format (v0.1)
		jsonData, err := json.Marshal(reads)
		if err != nil {
			return err
		}

		// Compress if enabled (for JSON only, binary writer handles its own compression)
		if w.compressor != nil && compressionType == "zstd" {
			compressed, err := w.compressor.Compress(jsonData)
			if err != nil {
				return fmt.Errorf("failed to compress chunk: %w", err)
			}
			data = compressed
		} else {
			data = jsonData
		}
	}

	if err := w.storage.WriteFile(chunkPath, data); err != nil {
		return err
	}

	// Calculate checksum (of compressed data)
	hash := sha256.Sum256(data)
	checksum := fmt.Sprintf("%x", hash)

	// Get file size
	fileSize := int64(len(data))

	// Record chunk info
	chunkInfo := ChunkInfo{
		Path:        filepath.Join("data", reference, filepath.Base(chunkPath)),
		Reference:   reference,
		Start:       start,
		End:         end,
		Reads:       len(reads),
		SizeBytes:   fileSize,
		Compression: compressionType,
		Checksum:    checksum,
		Created:     time.Now(),
	}

	if reference == "unmapped" {
		chunkInfo.Path = filepath.Join("data", "unmapped.chunk")
	}

	w.chunkInfos = append(w.chunkInfos, chunkInfo)

	return nil
}

// writeHeader writes the header file
func (w *Writer) writeHeader() error {
	data, err := json.MarshalIndent(w.header, "", "  ")
	if err != nil {
		return err
	}

	return w.storage.WriteFile("_header.json", data)
}

// writeMetadata writes the metadata file
func (w *Writer) writeMetadata() error {
	data, err := json.MarshalIndent(w.metadata, "", "  ")
	if err != nil {
		return err
	}

	return w.storage.WriteFile("_metadata.json", data)
}

// writeSpatialIndex writes the spatial index
func (w *Writer) writeSpatialIndex() error {
	// Build index from chunks
	index := SpatialIndex{
		References: make(map[string]ReferenceIndex),
	}

	for _, chunkInfo := range w.chunkInfos {
		if chunkInfo.Reference == "unmapped" {
			continue
		}

		refIndex, ok := index.References[chunkInfo.Reference]
		if !ok {
			// Find reference length from header
			refLength := 0
			for _, sq := range w.header.SQ {
				if sq["SN"] == chunkInfo.Reference {
					fmt.Sscanf(sq["LN"], "%d", &refLength)
					break
				}
			}

			refIndex = ReferenceIndex{
				Name:   chunkInfo.Reference,
				Length: refLength,
				Chunks: []ChunkInfo{},
			}
		}

		refIndex.Chunks = append(refIndex.Chunks, chunkInfo)
		index.References[chunkInfo.Reference] = refIndex
	}

	// Write index
	data, err := json.MarshalIndent(index, "", "  ")
	if err != nil {
		return err
	}

	return w.storage.WriteFile(filepath.Join("_index", "spatial.json"), data)
}
