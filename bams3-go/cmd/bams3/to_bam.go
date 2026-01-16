package main

import (
	"github.com/scttfrdmn/bams3-go/pkg/bam"
	"github.com/spf13/cobra"
)

var (
	toBamRegion string
	toBamNoIndex bool
)

var toBamCmd = &cobra.Command{
	Use:   "to-bam <input.bams3> <output.bam>",
	Short: "Convert BAMS3 dataset back to BAM format",
	Long: `Convert a BAMS3 dataset back to standard BAM format.

This enables compatibility with existing tools that require BAM files.
The conversion reconstructs a standard BAM file from BAMS3 chunks.

Examples:
  # Convert entire dataset to BAM
  bams3 to-bam sample.bams3 reconstructed.bam

  # Stream to stdout (for piping to other tools)
  bams3 to-bam sample.bams3 - | samtools view -c

  # Extract specific region
  bams3 to-bam sample.bams3 output.bam --region chr1:1000000-2000000

  # Stream region to GATK HaplotypeCaller
  bams3 to-bam sample.bams3 - --region chr17:41196312-41277500 | \
    gatk HaplotypeCaller -I /dev/stdin -R ref.fa -O output.vcf

Note: The output BAM will be coordinate-sorted if the BAMS3 dataset
was created from a coordinate-sorted BAM.`,
	Args: cobra.ExactArgs(2),
	RunE: func(cmd *cobra.Command, args []string) error {
		inputBAMS3 := args[0]
		outputBAM := args[1]

		opts := bam.ToBamOptions{
			Region:      toBamRegion,
			CreateIndex: !toBamNoIndex,
		}

		return bam.ConvertBAMS3ToBAM(inputBAMS3, outputBAM, opts)
	},
}

func init() {
	rootCmd.AddCommand(toBamCmd)

	toBamCmd.Flags().StringVarP(&toBamRegion, "region", "r", "",
		"Extract specific region (format: chr:start-end or chr)")
	toBamCmd.Flags().BoolVar(&toBamNoIndex, "no-index", false,
		"Skip creating BAM index file")
}
