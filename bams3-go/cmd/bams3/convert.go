package main

import (
	"github.com/scttfrdmn/bams3-go/pkg/bam"
	"github.com/spf13/cobra"
)

var (
	chunkSize int
)

var convertCmd = &cobra.Command{
	Use:   "convert <input.bam> <output.bams3>",
	Short: "Convert BAM file to BAMS3 format",
	Long: `Convert a standard BAM file to cloud-native BAMS3 format.

The BAMS3 format splits the BAM file into independent chunks (default 1Mbp)
that can be queried individually from S3 without downloading the entire file.

Example:
  bams3 convert sample.bam sample.bams3
  bams3 convert --chunk-size 5000000 large.bam large.bams3`,
	Args: cobra.ExactArgs(2),
	RunE: func(cmd *cobra.Command, args []string) error {
		inputBAM := args[0]
		outputBAMS3 := args[1]

		return bam.ConvertBAMToBAMS3(inputBAM, outputBAMS3, chunkSize)
	},
}

func init() {
	convertCmd.Flags().IntVar(&chunkSize, "chunk-size", 1000000,
		"Chunk size in base pairs (default 1Mbp)")
}
