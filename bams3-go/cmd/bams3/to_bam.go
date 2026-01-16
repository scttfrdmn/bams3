package main

import (
	"github.com/scttfrdmn/bams3-go/pkg/bam"
	"github.com/spf13/cobra"
)

var toBamCmd = &cobra.Command{
	Use:   "to-bam <input.bams3> <output.bam>",
	Short: "Convert BAMS3 dataset back to BAM format",
	Long: `Convert a BAMS3 dataset back to standard BAM format.

This enables compatibility with existing tools that require BAM files.
The conversion reconstructs a standard BAM file from BAMS3 chunks.

Example:
  bams3 to-bam sample.bams3 reconstructed.bam

Note: The output BAM will be coordinate-sorted if the BAMS3 dataset
was created from a coordinate-sorted BAM.`,
	Args: cobra.ExactArgs(2),
	RunE: func(cmd *cobra.Command, args []string) error {
		inputBAMS3 := args[0]
		outputBAM := args[1]

		return bam.ConvertBAMS3ToBAM(inputBAMS3, outputBAM)
	},
}

func init() {
	rootCmd.AddCommand(toBamCmd)
}
