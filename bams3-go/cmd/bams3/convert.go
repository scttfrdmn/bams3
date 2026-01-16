package main

import (
	"fmt"

	"github.com/scttfrdmn/bams3-go/pkg/bam"
	"github.com/scttfrdmn/bams3-go/pkg/bams3"
	"github.com/spf13/cobra"
)

var (
	chunkSizeStr string
	compression  string
	format       string
)

var convertCmd = &cobra.Command{
	Use:   "convert <input.bam> <output.bams3>",
	Short: "Convert BAM file to BAMS3 format",
	Long: `Convert a standard BAM file to cloud-native BAMS3 format.

The BAMS3 format splits the BAM file into independent chunks that can be
queried individually from S3 without downloading the entire file.

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

Examples:
  # Default settings (binary format, zstd compression, 1M chunks)
  bams3 convert sample.bam sample.bams3

  # Smaller chunks for low-latency queries
  bams3 convert --chunk-size 512K sample.bam sample.bams3

  # Legacy JSON format
  bams3 convert --format json sample.bam sample.bams3`,
	Args: cobra.ExactArgs(2),
	RunE: func(cmd *cobra.Command, args []string) error {
		inputBAM := args[0]
		outputBAMS3 := args[1]

		// Parse chunk size
		chunkSize, err := bams3.ParseChunkSize(chunkSizeStr)
		if err != nil {
			return fmt.Errorf("invalid chunk size: %w", err)
		}

		return bam.ConvertBAMToBAMS3(inputBAM, outputBAMS3, chunkSize, compression, format)
	},
}

func init() {
	convertCmd.Flags().StringVar(&chunkSizeStr, "chunk-size", "1M",
		"Chunk size: 256K, 512K, 1M, 2M, 4M, 8M (power of 2)")
	convertCmd.Flags().StringVar(&compression, "compression", "zstd",
		"Compression algorithm: none, zstd (default: zstd)")
	convertCmd.Flags().StringVar(&format, "format", "binary",
		"Chunk format: binary (v0.2), json (v0.1, legacy)")
}
