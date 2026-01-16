package bams3

import (
	"crypto/sha256"
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"time"
)

// Writer writes BAMS3 datasets
type Writer struct {
	path        string
	chunkSize   int
	compression string
	metadata    Metadata
	header      Header
	chunks      map[string][]ChunkRead // key: "ref:start:end"
	chunkInfos  []ChunkInfo
	statistics  Statistics
	compressor  *Compressor
}

// NewWriter creates a new BAMS3 writer
func NewWriter(path string, chunkSize int, compression string) (*Writer, error) {
	w := &Writer{
		path:        path,
		chunkSize:   chunkSize,
		compression: compression,
		chunks:      make(map[string][]ChunkRead),
		metadata: Metadata{
			Format:     "bams3",
			Version:    "0.1.0",
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

	// Create directory structure
	if err := os.MkdirAll(path, 0755); err != nil {
		return nil, fmt.Errorf("failed to create dataset directory: %w", err)
	}

	if err := os.MkdirAll(filepath.Join(path, "_index"), 0755); err != nil {
		return nil, fmt.Errorf("failed to create index directory: %w", err)
	}

	if err := os.MkdirAll(filepath.Join(path, "data"), 0755); err != nil {
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
		Name:  read.Name,
		Flag:  read.Flag,
		Ref:   read.ReferenceID,
		Pos:   read.Position,
		MapQ:  read.MappingQuality,
		CIGAR: read.CIGAR,
		Seq:   read.Sequence,
		Qual:  read.Quality,
		Tags:  read.Tags,
	}

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
		chunkPath = filepath.Join(w.path, "data", "unmapped.chunk")
	} else {
		refDir := filepath.Join(w.path, "data", reference)
		if err := os.MkdirAll(refDir, 0755); err != nil {
			return err
		}
		chunkPath = filepath.Join(refDir, fmt.Sprintf("%09d-%09d.chunk", start, end))
	}

	// Write chunk data
	data, err := json.Marshal(reads)
	if err != nil {
		return err
	}

	// Compress if enabled
	compressionType := w.compression
	if w.compressor != nil && compressionType == "zstd" {
		compressed, err := w.compressor.Compress(data)
		if err != nil {
			return fmt.Errorf("failed to compress chunk: %w", err)
		}
		data = compressed
	}

	if err := os.WriteFile(chunkPath, data, 0644); err != nil {
		return err
	}

	// Calculate checksum (of compressed data)
	hash := sha256.Sum256(data)
	checksum := fmt.Sprintf("%x", hash)

	// Get file size
	info, err := os.Stat(chunkPath)
	if err != nil {
		return err
	}

	// Record chunk info
	chunkInfo := ChunkInfo{
		Path:        filepath.Join("data", reference, filepath.Base(chunkPath)),
		Reference:   reference,
		Start:       start,
		End:         end,
		Reads:       len(reads),
		SizeBytes:   info.Size(),
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
	path := filepath.Join(w.path, "_header.json")
	data, err := json.MarshalIndent(w.header, "", "  ")
	if err != nil {
		return err
	}

	return os.WriteFile(path, data, 0644)
}

// writeMetadata writes the metadata file
func (w *Writer) writeMetadata() error {
	path := filepath.Join(w.path, "_metadata.json")
	data, err := json.MarshalIndent(w.metadata, "", "  ")
	if err != nil {
		return err
	}

	return os.WriteFile(path, data, 0644)
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
	path := filepath.Join(w.path, "_index", "spatial.json")
	data, err := json.MarshalIndent(index, "", "  ")
	if err != nil {
		return err
	}

	return os.WriteFile(path, data, 0644)
}
