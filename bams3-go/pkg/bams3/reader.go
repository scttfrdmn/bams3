package bams3

import (
	"encoding/json"
	"fmt"
	"path/filepath"
	"strings"
)

// Reader reads BAMS3 datasets
type Reader struct {
	dataset Dataset
	storage Storage // Storage backend (local or S3)
}

// OpenDataset opens a BAMS3 dataset
func OpenDataset(path string) (*Reader, error) {
	// Create storage backend (auto-detects local vs S3)
	storage, err := NewStorage(path)
	if err != nil {
		return nil, fmt.Errorf("failed to create storage backend: %w", err)
	}

	r := &Reader{
		storage: storage,
	}

	// Load metadata
	if err := r.loadMetadata("_metadata.json"); err != nil {
		return nil, fmt.Errorf("failed to load metadata: %w", err)
	}

	// Load header
	if err := r.loadHeader("_header.json"); err != nil {
		return nil, fmt.Errorf("failed to load header: %w", err)
	}

	// Load spatial index
	if err := r.loadIndex(filepath.Join("_index", "spatial.json")); err != nil {
		return nil, fmt.Errorf("failed to load index: %w", err)
	}

	r.dataset.Path = path

	return r, nil
}

// loadMetadata loads the metadata file
func (r *Reader) loadMetadata(path string) error {
	data, err := r.storage.ReadFile(path)
	if err != nil {
		return err
	}

	return json.Unmarshal(data, &r.dataset.Metadata)
}

// loadHeader loads the header file
func (r *Reader) loadHeader(path string) error {
	data, err := r.storage.ReadFile(path)
	if err != nil {
		return err
	}

	return json.Unmarshal(data, &r.dataset.Header)
}

// loadIndex loads the spatial index
func (r *Reader) loadIndex(path string) error {
	data, err := r.storage.ReadFile(path)
	if err != nil {
		return err
	}

	return json.Unmarshal(data, &r.dataset.Index)
}

// GetMetadata returns the dataset metadata
func (r *Reader) GetMetadata() Metadata {
	return r.dataset.Metadata
}

// GetHeader returns the SAM header
func (r *Reader) GetHeader() Header {
	return r.dataset.Header
}

// GetStatistics returns dataset statistics
func (r *Reader) GetStatistics() Statistics {
	return r.dataset.Metadata.Statistics
}

// FindChunks finds chunks that overlap with a region
func (r *Reader) FindChunks(region Region) ([]ChunkInfo, error) {
	var overlapping []ChunkInfo

	for _, chunk := range r.dataset.Metadata.Chunks {
		if chunk.Reference != region.Reference {
			continue
		}

		// Check for overlap
		if chunk.End > region.Start && chunk.Start < region.End {
			overlapping = append(overlapping, chunk)
		}
	}

	return overlapping, nil
}

// QueryRegion queries reads in a specific genomic region
func (r *Reader) QueryRegion(region Region) ([]Read, error) {
	// Find overlapping chunks
	chunks, err := r.FindChunks(region)
	if err != nil {
		return nil, err
	}

	var reads []Read

	// Load and filter reads from each chunk
	for _, chunkInfo := range chunks {
		chunkReads, err := r.loadChunk(chunkInfo.Path)
		if err != nil {
			return nil, fmt.Errorf("failed to load chunk %s: %w", chunkInfo.Path, err)
		}

		// Filter reads to exact region
		for _, read := range chunkReads {
			if read.Position >= region.Start && read.Position < region.End {
				reads = append(reads, read)
			}
		}
	}

	return reads, nil
}

// loadChunk loads a chunk file and converts to Read objects
func (r *Reader) loadChunk(path string) ([]Read, error) {
	data, err := r.storage.ReadFile(path)
	if err != nil {
		return nil, err
	}

	// Find chunk info to check compression
	var chunkInfo *ChunkInfo
	for i := range r.dataset.Metadata.Chunks {
		if r.dataset.Metadata.Chunks[i].Path == path {
			chunkInfo = &r.dataset.Metadata.Chunks[i]
			break
		}
	}

	// Decompress if needed (for both binary and JSON formats)
	if chunkInfo != nil && chunkInfo.Compression == "zstd" {
		decompressed, err := Decompress(data)
		if err != nil {
			return nil, fmt.Errorf("failed to decompress chunk: %w", err)
		}
		data = decompressed
	}

	// Auto-detect format by checking for binary magic number (after decompression)
	// Magic number 0x42414D33 is stored in little-endian as: 0x33, 0x4D, 0x41, 0x42
	isBinary := len(data) >= 4 &&
		data[0] == 0x33 && data[1] == 0x4D && data[2] == 0x41 && data[3] == 0x42

	var chunkReads []ChunkRead

	if isBinary {
		// Use binary reader
		binaryReader := &BinaryChunkReader{}
		var err error
		chunkReads, err = binaryReader.ReadChunk(data)
		if err != nil {
			return nil, fmt.Errorf("failed to read binary chunk: %w", err)
		}
	} else {
		// Use JSON reader
		if err := json.Unmarshal(data, &chunkReads); err != nil {
			return nil, fmt.Errorf("failed to parse JSON chunk: %w", err)
		}
	}

	// Convert ChunkRead to Read
	reads := make([]Read, len(chunkReads))
	for i, cr := range chunkReads {
		reads[i] = Read{
			Name:           cr.Name,
			Flag:           cr.Flag,
			ReferenceID:    cr.ReferenceID,
			Position:       cr.Position,
			MappingQuality: cr.MappingQuality,
			CIGAR:          cr.CIGAR,
			CIGAROps:       cr.CIGAROps,
			Sequence:       cr.Sequence,
			Quality:        cr.Quality,
			Tags:           cr.Tags,
		}

		// Fall back to legacy fields if new fields are empty
		if reads[i].ReferenceID == 0 && cr.Ref != 0 {
			reads[i].ReferenceID = cr.Ref
		}
		if reads[i].Position == 0 && cr.Pos != 0 {
			reads[i].Position = cr.Pos
		}
		if reads[i].MappingQuality == 0 && cr.MapQ != 0 {
			reads[i].MappingQuality = cr.MapQ
		}
		if reads[i].Sequence == "" && cr.Seq != "" {
			reads[i].Sequence = cr.Seq
		}
		if reads[i].Quality == "" && cr.Qual != "" {
			reads[i].Quality = cr.Qual
		}
	}

	return reads, nil
}

// ParseRegion parses a region string like "chr1:1000000-2000000"
func ParseRegion(regionStr string) (Region, error) {
	region := Region{}

	// Split on ':'
	parts := strings.Split(regionStr, ":")
	if len(parts) != 2 {
		return region, fmt.Errorf("invalid region format: %s (expected chr:start-end)", regionStr)
	}

	region.Reference = parts[0]

	// Split on '-'
	posParts := strings.Split(parts[1], "-")
	if len(posParts) != 2 {
		return region, fmt.Errorf("invalid region format: %s (expected chr:start-end)", regionStr)
	}

	_, err := fmt.Sscanf(posParts[0], "%d", &region.Start)
	if err != nil {
		return region, fmt.Errorf("invalid start position: %w", err)
	}

	_, err = fmt.Sscanf(posParts[1], "%d", &region.End)
	if err != nil {
		return region, fmt.Errorf("invalid end position: %w", err)
	}

	return region, nil
}

// Close closes the reader
func (r *Reader) Close() error {
	// Nothing to close for local files
	return nil
}
