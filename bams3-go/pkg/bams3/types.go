package bams3

import "time"

// Dataset represents a BAMS3 dataset
type Dataset struct {
	Path     string
	Metadata Metadata
	Header   Header
	Index    SpatialIndex
}

// Metadata contains dataset metadata
type Metadata struct {
	Format      string            `json:"format"`
	Version     string            `json:"version"`
	Created     time.Time         `json:"created"`
	CreatedBy   string            `json:"created_by"`
	Source      Source            `json:"source"`
	Statistics  Statistics        `json:"statistics"`
	Chunks      []ChunkInfo       `json:"chunks"`
	ChunkSize   int               `json:"chunk_size"`
	Compression CompressionConfig `json:"compression"`
}

// Source describes the original data source
type Source struct {
	File   string `json:"file"`
	Format string `json:"format"`
}

// Statistics contains dataset statistics
type Statistics struct {
	TotalReads     int     `json:"total_reads"`
	MappedReads    int     `json:"mapped_reads"`
	UnmappedReads  int     `json:"unmapped_reads"`
	DuplicateReads int     `json:"duplicate_reads"`
	TotalBases     int64   `json:"total_bases"`
	MeanCoverage   float64 `json:"mean_coverage"`
}

// ChunkInfo describes a chunk
type ChunkInfo struct {
	Path        string    `json:"path"`
	Reference   string    `json:"reference"`
	Start       int       `json:"start"`
	End         int       `json:"end"`
	Reads       int       `json:"reads"`
	SizeBytes   int64     `json:"size_bytes"`
	Compression string    `json:"compression"`
	Checksum    string    `json:"checksum"`
	Created     time.Time `json:"created"`
}

// CompressionConfig describes compression settings
type CompressionConfig struct {
	Algorithm string `json:"algorithm"`
}

// Header represents SAM header
type Header struct {
	HD map[string]string   `json:"HD"` // Header line
	SQ []map[string]string `json:"SQ"` // Reference sequences
	RG []map[string]string `json:"RG"` // Read groups
	PG []map[string]string `json:"PG"` // Programs
	CO []string            `json:"CO"` // Comments
}

// SpatialIndex maps genomic regions to chunks
type SpatialIndex struct {
	References map[string]ReferenceIndex `json:"references"`
}

// ReferenceIndex contains chunks for a reference
type ReferenceIndex struct {
	Name   string      `json:"name"`
	Length int         `json:"length"`
	Chunks []ChunkInfo `json:"chunks"`
}

// Read represents an alignment read
type Read struct {
	Name            string
	Flag            int
	ReferenceID     int
	Position        int
	MappingQuality  uint8
	CIGAR           string
	CIGAROps        []CIGAROperation // For binary format
	MateReferenceID int
	MatePosition    int
	TemplateLength  int
	Sequence        string
	Quality         string
	Tags            map[string]interface{}
}

// CIGAROperation represents a single CIGAR operation
type CIGAROperation struct {
	Type   byte // MIDNSHP=X
	Length int
}

// ChunkData represents the contents of a chunk
type ChunkData struct {
	Reference string      `json:"reference"`
	Start     int         `json:"start"`
	End       int         `json:"end"`
	Reads     []ChunkRead `json:"reads"`
}

// ChunkRead is the representation of a read in a chunk (supports both JSON and binary)
type ChunkRead struct {
	// Common fields
	Name           string                 `json:"name"`
	Flag           int                    `json:"flag"`
	ReferenceID    int                    `json:"reference_id,omitempty"`    // Binary format
	Position       int                    `json:"position,omitempty"`        // Binary format
	MappingQuality uint8                  `json:"mapping_quality,omitempty"` // Binary format
	CIGAR          string                 `json:"cigar"`
	CIGAROps       []CIGAROperation       `json:"cigar_ops,omitempty"` // Binary format
	Sequence       string                 `json:"sequence,omitempty"`  // Binary format
	Quality        string                 `json:"quality,omitempty"`   // Binary format
	Tags           map[string]interface{} `json:"tags,omitempty"`

	// Legacy JSON fields (v0.1)
	Ref  int    `json:"ref,omitempty"`
	Pos  int    `json:"pos,omitempty"`
	MapQ uint8  `json:"mapq,omitempty"`
	Seq  string `json:"seq,omitempty"`
	Qual string `json:"qual,omitempty"`
}

// Region represents a genomic region query
type Region struct {
	Reference string
	Start     int
	End       int
}
