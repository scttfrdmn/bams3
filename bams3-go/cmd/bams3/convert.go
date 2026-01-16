package main

import (
	"fmt"
	"os"

	"github.com/scttfrdmn/bams3-go/pkg/bam"
	"github.com/scttfrdmn/bams3-go/pkg/bams3"
	"github.com/spf13/cobra"
)

var (
	chunkSizeStr string
	compression  string
	format       string
	workers      int
	fromStdin    bool
	showConfig   bool
	bufferSize   string
	sortBuffer   string
)

var convertCmd = &cobra.Command{
	Use:   "convert [input.bam] <output.bams3>",
	Short: "Convert BAM file to BAMS3 format",
	Long: `Convert a standard BAM file to cloud-native BAMS3 format.

The BAMS3 format splits the BAM file into independent chunks that can be
queried individually from S3 without downloading the entire file.

Streaming Conversion:
  Use --stdin to read SAM/BAM from stdin, enabling zero-copy pipelines:
    bwa mem ref.fa reads.fq | bams3 convert --stdin output.bams3
    samtools view input.bam | bams3 convert --stdin output.bams3

  This eliminates intermediate files and enables direct upload to S3:
    bwa mem ref.fa reads.fq | bams3 convert --stdin s3://bucket/output.bams3

Chunk Size:
  Power-of-2 sizes optimize S3 access patterns:
    256K - Ultra-low latency queries
    512K - Low latency queries
    1M   - Default, balanced performance (recommended)
    2M   - Larger queries, moderate access
    4M   - Archival, infrequent access
    8M   - Cold storage

Format:
  binary - Compact binary format (v0.2.0, default)
  json   - Human-readable JSON format (v0.1.0, legacy)

Smart Defaults:
  Workers: Auto-detected from CPU count
  Buffer: 25% of available RAM
  Sort buffer: min(8GB, 25% of RAM)
  All settings can be overridden with flags

Examples:
  # File conversion (traditional)
  bams3 convert sample.bam sample.bams3

  # Streaming conversion from aligner
  bwa mem ref.fa reads.fq | bams3 convert --stdin output.bams3

  # Direct S3 upload while converting
  bwa mem ref.fa reads.fq | bams3 convert --stdin s3://bucket/output.bams3

  # Custom settings
  bwa mem ref.fa reads.fq | bams3 convert --stdin output.bams3 \
    --workers 16 --buffer 8G --sort-buffer 16G

  # Show effective configuration
  bams3 convert --stdin --show-config output.bams3 < input.sam`,
	Args: cobra.RangeArgs(1, 2),
	RunE: func(cmd *cobra.Command, args []string) error {
		// Parse chunk size
		chunkSize, err := bams3.ParseChunkSize(chunkSizeStr)
		if err != nil {
			return fmt.Errorf("invalid chunk size: %w", err)
		}

		// Handle --stdin mode
		if fromStdin {
			if len(args) != 1 {
				return fmt.Errorf("--stdin mode requires exactly one argument (output path)")
			}
			return convertFromStdin(args[0], chunkSize)
		}

		// Handle file mode
		if len(args) != 2 {
			return fmt.Errorf("file mode requires two arguments (input.bam output.bams3)")
		}

		inputBAM := args[0]
		outputBAMS3 := args[1]

		return bam.ConvertBAMToBAMS3(inputBAM, outputBAMS3, chunkSize, compression, format, workers)
	},
}

func init() {
	convertCmd.Flags().StringVar(&chunkSizeStr, "chunk-size", "1M",
		"Chunk size: 256K, 512K, 1M, 2M, 4M, 8M (power of 2)")
	convertCmd.Flags().StringVar(&compression, "compression", "zstd",
		"Compression algorithm: none, zstd (default: zstd)")
	convertCmd.Flags().StringVar(&format, "format", "binary",
		"Chunk format: binary (v0.2), json (v0.1, legacy)")
	convertCmd.Flags().IntVar(&workers, "workers", 0,
		"Number of parallel workers (0 = auto-detect CPU count, 1 = sequential)")

	// Streaming mode flags
	convertCmd.Flags().BoolVar(&fromStdin, "stdin", false,
		"Read SAM/BAM from stdin (enables zero-copy pipelines)")
	convertCmd.Flags().BoolVar(&showConfig, "show-config", false,
		"Show effective configuration (workers, buffers, etc.)")
	convertCmd.Flags().StringVar(&bufferSize, "buffer", "",
		"Read buffer size (e.g., 4G, 8G) - default: 25% of RAM")
	convertCmd.Flags().StringVar(&sortBuffer, "sort-buffer", "",
		"Sort buffer size (e.g., 8G, 16G) - default: min(8G, 25% of RAM)")
}

// convertFromStdin handles streaming conversion from stdin
func convertFromStdin(outputPath string, chunkSize int) error {
	// Create stream config with defaults
	config := bams3.NewStreamConfig()
	config.ChunkSize = chunkSize
	config.Compression = compression
	config.Format = format

	// Override workers if specified
	if workers > 0 {
		config.Workers = workers
	}

	// Parse buffer sizes if specified
	if bufferSize != "" {
		size, err := bams3.ParseSize(bufferSize)
		if err != nil {
			return fmt.Errorf("invalid buffer size: %w", err)
		}
		config.BufferSize = size
	}

	if sortBuffer != "" {
		size, err := bams3.ParseSize(sortBuffer)
		if err != nil {
			return fmt.Errorf("invalid sort buffer size: %w", err)
		}
		config.SortBufferSize = size
	}

	// Show config if requested
	if showConfig {
		config.ShowConfig()
		return nil
	}

	// Create converter
	converter, err := bams3.NewStreamConverter(outputPath, config)
	if err != nil {
		return fmt.Errorf("failed to create converter: %w", err)
	}

	// Source info
	source := bams3.Source{
		File:   "stdin",
		Format: "SAM/BAM",
	}

	// Convert from stdin (header will be parsed from input)
	return converter.ConvertFromStream(os.Stdin, source)
}
