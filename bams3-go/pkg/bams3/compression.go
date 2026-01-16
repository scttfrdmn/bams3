package bams3

import (
	"bytes"
	"fmt"
	"io"

	"github.com/klauspost/compress/zstd"
)

// Compressor handles data compression
type Compressor struct {
	encoder *zstd.Encoder
	decoder *zstd.Decoder
}

// NewCompressor creates a new compressor
func NewCompressor() (*Compressor, error) {
	encoder, err := zstd.NewWriter(nil)
	if err != nil {
		return nil, fmt.Errorf("failed to create zstd encoder: %w", err)
	}

	decoder, err := zstd.NewReader(nil)
	if err != nil {
		return nil, fmt.Errorf("failed to create zstd decoder: %w", err)
	}

	return &Compressor{
		encoder: encoder,
		decoder: decoder,
	}, nil
}

// Compress compresses data using zstd
func (c *Compressor) Compress(data []byte) ([]byte, error) {
	return c.encoder.EncodeAll(data, make([]byte, 0, len(data))), nil
}

// Decompress decompresses data using zstd
func (c *Compressor) Decompress(data []byte) ([]byte, error) {
	return c.decoder.DecodeAll(data, nil)
}

// CompressWithLevel compresses data with a specific level
func CompressWithLevel(data []byte, level int) ([]byte, error) {
	var encoderLevel zstd.EncoderLevel
	switch level {
	case 1:
		encoderLevel = zstd.SpeedFastest
	case 2:
		encoderLevel = zstd.SpeedDefault
	case 3:
		encoderLevel = zstd.SpeedBetterCompression
	default:
		encoderLevel = zstd.SpeedDefault
	}

	encoder, err := zstd.NewWriter(nil, zstd.WithEncoderLevel(encoderLevel))
	if err != nil {
		return nil, fmt.Errorf("failed to create encoder: %w", err)
	}

	return encoder.EncodeAll(data, make([]byte, 0, len(data))), nil
}

// Decompress decompresses zstd-compressed data
func Decompress(data []byte) ([]byte, error) {
	decoder, err := zstd.NewReader(bytes.NewReader(data))
	if err != nil {
		return nil, fmt.Errorf("failed to create decoder: %w", err)
	}
	defer decoder.Close()

	return io.ReadAll(decoder)
}

// Close closes the compressor
func (c *Compressor) Close() error {
	if c.encoder != nil {
		c.encoder.Close()
	}
	if c.decoder != nil {
		c.decoder.Close()
	}
	return nil
}
