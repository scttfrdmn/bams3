package bams3

import (
	"encoding/binary"
	"fmt"
	"io"
	"math"
	"strconv"
	"strings"
)

// Binary format constants
const (
	BinaryMagic   uint32 = 0x42414D33 // "BAM3"
	BinaryVersion uint16 = 0x0200     // v0.2.0
)

// Compression flags
const (
	CompressionNone uint16 = 0x0
	CompressionZstd uint16 = 0x1
)

// Base encoding (4-bit per base)
const (
	BaseA byte = 1
	BaseC byte = 2
	BaseG byte = 4
	BaseT byte = 8
	BaseN byte = 15
)

// BinaryChunkWriter writes reads in binary format
type BinaryChunkWriter struct {
	compression string
	compressor  *Compressor
}

// NewBinaryChunkWriter creates a new binary chunk writer
func NewBinaryChunkWriter(compression string) (*BinaryChunkWriter, error) {
	w := &BinaryChunkWriter{
		compression: compression,
	}

	// Initialize compressor if needed
	if compression == "zstd" {
		compressor, err := NewCompressor()
		if err != nil {
			return nil, fmt.Errorf("failed to create compressor: %w", err)
		}
		w.compressor = compressor
	}

	return w, nil
}

// WriteChunk writes a chunk of reads in binary format
func (w *BinaryChunkWriter) WriteChunk(reads []ChunkRead) ([]byte, error) {
	// Estimate buffer size (average ~150 bytes per read)
	estimatedSize := 16 + len(reads)*150
	buf := make([]byte, 0, estimatedSize)

	// Write chunk header
	buf = w.writeChunkHeader(buf, len(reads))

	// Write each record
	for _, read := range reads {
		recordBuf, err := w.encodeRecord(read)
		if err != nil {
			return nil, fmt.Errorf("failed to encode record: %w", err)
		}
		buf = append(buf, recordBuf...)
	}

	// Compress if enabled
	if w.compressor != nil && w.compression == "zstd" {
		compressed, err := w.compressor.Compress(buf)
		if err != nil {
			return nil, fmt.Errorf("failed to compress chunk: %w", err)
		}
		return compressed, nil
	}

	return buf, nil
}

// writeChunkHeader writes the 16-byte chunk header
func (w *BinaryChunkWriter) writeChunkHeader(buf []byte, numRecords int) []byte {
	header := make([]byte, 16)

	// Magic number
	binary.LittleEndian.PutUint32(header[0:4], BinaryMagic)

	// Version
	binary.LittleEndian.PutUint16(header[4:6], BinaryVersion)

	// Flags
	var flags uint16 = CompressionNone
	if w.compression == "zstd" {
		flags = CompressionZstd
	}
	binary.LittleEndian.PutUint16(header[6:8], flags)

	// Number of records
	binary.LittleEndian.PutUint32(header[8:12], uint32(numRecords))

	// Reserved
	binary.LittleEndian.PutUint32(header[12:16], 0)

	return append(buf, header...)
}

// encodeRecord encodes a single read record
func (w *BinaryChunkWriter) encodeRecord(read ChunkRead) ([]byte, error) {
	// Calculate sizes
	nameLen := len(read.Name)
	seqLen := len(read.Sequence)
	cigarOps := len(read.CIGAROps)
	numTags := len(read.Tags)

	// Estimate record size
	recordSize := 12 + // header
		nameLen + // name
		4 + // reference_id
		4 + // position
		1 + // mapping_quality
		2 + // flags
		cigarOps*4 + // cigar
		(seqLen+1)/2 + // sequence (2 bases per byte)
		seqLen // quality

	// Add tag sizes (estimate)
	for _, tagValue := range read.Tags {
		recordSize += 2 + 1 + len(fmt.Sprintf("%v", tagValue))
	}

	buf := make([]byte, 0, recordSize)

	// Record header (12 bytes) - we'll update record_size at the end
	headerStart := len(buf)
	buf = append(buf, make([]byte, 12)...)

	// Store sizes in header
	binary.LittleEndian.PutUint16(buf[headerStart+4:], uint16(nameLen))
	binary.LittleEndian.PutUint16(buf[headerStart+6:], uint16(seqLen))
	binary.LittleEndian.PutUint16(buf[headerStart+8:], uint16(cigarOps))
	binary.LittleEndian.PutUint16(buf[headerStart+10:], uint16(numTags))

	// Read name
	buf = append(buf, []byte(read.Name)...)

	// Reference ID
	refID := int32(read.ReferenceID)
	binary.LittleEndian.PutUint32(buf[len(buf):len(buf)+4], uint32(refID))
	buf = buf[:len(buf)+4]

	// Position
	pos := int32(read.Position)
	binary.LittleEndian.PutUint32(buf[len(buf):len(buf)+4], uint32(pos))
	buf = buf[:len(buf)+4]

	// Mapping quality
	buf = append(buf, read.MappingQuality)

	// Flags
	flags := uint16(read.Flag)
	binary.LittleEndian.PutUint16(buf[len(buf):len(buf)+2], flags)
	buf = buf[:len(buf)+2]

	// CIGAR operations
	for _, op := range read.CIGAROps {
		// Pack: upper 28 bits = length, lower 4 bits = op type
		packed := (uint32(op.Length) << 4) | uint32(op.Type)
		binary.LittleEndian.PutUint32(buf[len(buf):len(buf)+4], packed)
		buf = buf[:len(buf)+4]
	}

	// Sequence (4-bit encoding, 2 bases per byte)
	seqBytes := encodeSequence(read.Sequence)
	buf = append(buf, seqBytes...)

	// Quality scores
	buf = append(buf, []byte(read.Quality)...)

	// Tags (simplified for now - store as strings)
	for tagName, tagValue := range read.Tags {
		// Tag name (2 bytes)
		if len(tagName) >= 2 {
			buf = append(buf, tagName[0], tagName[1])
		} else {
			buf = append(buf, tagName[0], 0)
		}

		// Type (1 byte) - 'Z' for string
		buf = append(buf, 'Z')

		// Value as string
		valueStr := fmt.Sprintf("%v", tagValue)
		buf = append(buf, []byte(valueStr)...)
		buf = append(buf, 0) // null terminator
	}

	// Update record size in header
	actualSize := len(buf)
	binary.LittleEndian.PutUint32(buf[headerStart:headerStart+4], uint32(actualSize))

	return buf, nil
}

// encodeSequence encodes DNA sequence to 4-bit format (2 bases per byte)
func encodeSequence(seq string) []byte {
	numBytes := (len(seq) + 1) / 2
	encoded := make([]byte, numBytes)

	for i := 0; i < len(seq); i++ {
		base := encodeBase(seq[i])
		byteIdx := i / 2
		if i%2 == 0 {
			// Store in upper 4 bits
			encoded[byteIdx] = base << 4
		} else {
			// Store in lower 4 bits
			encoded[byteIdx] |= base
		}
	}

	return encoded
}

// encodeBase converts a base character to 4-bit encoding
func encodeBase(b byte) byte {
	switch b {
	case 'A', 'a':
		return BaseA
	case 'C', 'c':
		return BaseC
	case 'G', 'g':
		return BaseG
	case 'T', 't':
		return BaseT
	case 'N', 'n':
		return BaseN
	default:
		return 0
	}
}

// ParseSize parses size string (e.g., "1M", "512K", "8G") to bytes
func ParseSize(sizeStr string) (int64, error) {
	sizeStr = strings.ToUpper(strings.TrimSpace(sizeStr))

	var multiplier int64 = 1
	if strings.HasSuffix(sizeStr, "K") {
		multiplier = 1024
		sizeStr = sizeStr[:len(sizeStr)-1]
	} else if strings.HasSuffix(sizeStr, "M") {
		multiplier = 1024 * 1024
		sizeStr = sizeStr[:len(sizeStr)-1]
	} else if strings.HasSuffix(sizeStr, "G") {
		multiplier = 1024 * 1024 * 1024
		sizeStr = sizeStr[:len(sizeStr)-1]
	}

	value, err := strconv.ParseInt(sizeStr, 10, 64)
	if err != nil {
		return 0, fmt.Errorf("invalid size: %s", sizeStr)
	}

	return value * multiplier, nil
}

// ParseChunkSize parses chunk size string (e.g., "1M", "512K") to bytes
func ParseChunkSize(sizeStr string) (int, error) {
	if len(sizeStr) < 2 {
		return 0, fmt.Errorf("invalid chunk size format: %s", sizeStr)
	}

	// Extract number and unit
	var value float64
	var unit string

	_, err := fmt.Sscanf(sizeStr, "%f%s", &value, &unit)
	if err != nil {
		return 0, fmt.Errorf("failed to parse chunk size: %w", err)
	}

	// Convert to bytes
	var multiplier int64
	switch unit {
	case "K", "KB":
		multiplier = 1024
	case "M", "MB":
		multiplier = 1024 * 1024
	case "G", "GB":
		multiplier = 1024 * 1024 * 1024
	default:
		return 0, fmt.Errorf("unknown size unit: %s", unit)
	}

	bytes := int64(value * float64(multiplier))

	// Validate it's a power of 2
	if !isPowerOfTwo(bytes) {
		return 0, fmt.Errorf("chunk size must be a power of 2, got: %d", bytes)
	}

	// Validate range (256KB to 8MB)
	if bytes < 256*1024 || bytes > 8*1024*1024 {
		return 0, fmt.Errorf("chunk size must be between 256K and 8M, got: %d", bytes)
	}

	return int(bytes), nil
}

// isPowerOfTwo checks if a number is a power of 2
func isPowerOfTwo(n int64) bool {
	return n > 0 && (n&(n-1)) == 0
}

// FormatChunkSize formats bytes as human-readable size
func FormatChunkSize(bytes int) string {
	if bytes >= 1024*1024 {
		mb := float64(bytes) / (1024 * 1024)
		if mb == math.Floor(mb) {
			return fmt.Sprintf("%dM", int(mb))
		}
		return fmt.Sprintf("%.1fM", mb)
	} else if bytes >= 1024 {
		kb := float64(bytes) / 1024
		if kb == math.Floor(kb) {
			return fmt.Sprintf("%dK", int(kb))
		}
		return fmt.Sprintf("%.1fK", kb)
	}
	return fmt.Sprintf("%d", bytes)
}

// BinaryChunkReader reads chunks in binary format
type BinaryChunkReader struct{}

// ReadChunk reads a chunk in binary format
func (r *BinaryChunkReader) ReadChunk(data []byte) ([]ChunkRead, error) {
	if len(data) < 16 {
		return nil, fmt.Errorf("chunk too small: %d bytes", len(data))
	}

	// Read chunk header
	magic := binary.LittleEndian.Uint32(data[0:4])
	if magic != BinaryMagic {
		return nil, fmt.Errorf("invalid magic number: 0x%08X", magic)
	}

	version := binary.LittleEndian.Uint16(data[4:6])
	if version != BinaryVersion {
		return nil, fmt.Errorf("unsupported version: 0x%04X", version)
	}

	// flags := binary.LittleEndian.Uint16(data[6:8])  // flags not used currently
	numRecords := binary.LittleEndian.Uint32(data[8:12])

	// Note: Decompression is handled by the caller (reader.go)
	// The data passed here should already be decompressed

	// Read records
	reads := make([]ChunkRead, 0, numRecords)
	offset := 16 // Skip header

	for i := uint32(0); i < numRecords; i++ {
		read, bytesRead, err := r.decodeRecord(data[offset:])
		if err != nil {
			return nil, fmt.Errorf("failed to decode record %d: %w", i, err)
		}
		reads = append(reads, read)
		offset += bytesRead
	}

	return reads, nil
}

// decodeRecord decodes a single record
func (r *BinaryChunkReader) decodeRecord(data []byte) (ChunkRead, int, error) {
	if len(data) < 12 {
		return ChunkRead{}, 0, fmt.Errorf("record too small")
	}

	// Read record header
	recordSize := binary.LittleEndian.Uint32(data[0:4])
	nameLen := binary.LittleEndian.Uint16(data[4:6])
	seqLen := binary.LittleEndian.Uint16(data[6:8])
	cigarOps := binary.LittleEndian.Uint16(data[8:10])
	numTags := binary.LittleEndian.Uint16(data[10:12])

	offset := 12
	read := ChunkRead{}

	// Read name
	read.Name = string(data[offset : offset+int(nameLen)])
	offset += int(nameLen)

	// Reference ID
	read.ReferenceID = int(int32(binary.LittleEndian.Uint32(data[offset : offset+4])))
	offset += 4

	// Position
	read.Position = int(int32(binary.LittleEndian.Uint32(data[offset : offset+4])))
	offset += 4

	// Mapping quality
	read.MappingQuality = data[offset]
	offset++

	// Flags
	read.Flag = int(binary.LittleEndian.Uint16(data[offset : offset+2]))
	offset += 2

	// CIGAR operations
	read.CIGAROps = make([]CIGAROperation, cigarOps)
	for i := 0; i < int(cigarOps); i++ {
		packed := binary.LittleEndian.Uint32(data[offset : offset+4])
		read.CIGAROps[i] = CIGAROperation{
			Length: int(packed >> 4),
			Type:   byte(packed & 0xF),
		}
		offset += 4
	}

	// Sequence
	seqBytes := (int(seqLen) + 1) / 2
	read.Sequence = decodeSequence(data[offset:offset+seqBytes], int(seqLen))
	offset += seqBytes

	// Quality
	read.Quality = string(data[offset : offset+int(seqLen)])
	offset += int(seqLen)

	// Tags
	read.Tags = make(map[string]interface{})
	for i := 0; i < int(numTags); i++ {
		// Tag name (2 bytes)
		tagName := string(data[offset : offset+2])
		offset += 2

		// Type (1 byte)
		tagType := data[offset]
		offset++

		// Value (null-terminated string for now)
		if tagType == 'Z' {
			endIdx := offset
			for endIdx < len(data) && data[endIdx] != 0 {
				endIdx++
			}
			read.Tags[tagName] = string(data[offset:endIdx])
			offset = endIdx + 1
		}
	}

	return read, int(recordSize), nil
}

// decodeSequence decodes 4-bit encoded sequence
func decodeSequence(data []byte, length int) string {
	bases := make([]byte, length)
	for i := 0; i < length; i++ {
		byteIdx := i / 2
		var encoded byte
		if i%2 == 0 {
			// Upper 4 bits
			encoded = data[byteIdx] >> 4
		} else {
			// Lower 4 bits
			encoded = data[byteIdx] & 0xF
		}
		bases[i] = decodeBase(encoded)
	}
	return string(bases)
}

// decodeBase converts 4-bit encoding to base character
func decodeBase(b byte) byte {
	switch b {
	case BaseA:
		return 'A'
	case BaseC:
		return 'C'
	case BaseG:
		return 'G'
	case BaseT:
		return 'T'
	case BaseN:
		return 'N'
	default:
		return 'N'
	}
}

// Ensure io.Writer interface (not actually used as writer, but for consistency)
var _ io.Writer = (*BinaryChunkWriter)(nil)

func (w *BinaryChunkWriter) Write(p []byte) (n int, err error) {
	// Not used - WriteChunk is the main method
	return 0, fmt.Errorf("use WriteChunk method instead")
}
