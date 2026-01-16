package main

import (
	"fmt"

	"github.com/scttfrdmn/bams3-go/pkg/bams3"
	"github.com/spf13/cobra"
)

var (
	countOnly bool
	showReads int
)

var queryCmd = &cobra.Command{
	Use:   "query <dataset.bams3> <region>",
	Short: "Query reads from a BAMS3 dataset",
	Long: `Query reads from a specific genomic region in a BAMS3 dataset.

The region format is: chr:start-end (e.g., chr1:1000000-2000000)

Only the chunks overlapping the query region are loaded, making this
much faster than downloading and querying the entire BAM file.

Examples:
  bams3 query sample.bams3 chr1:1000000-2000000
  bams3 query sample.bams3 chr1:1000000-2000000 --count
  bams3 query sample.bams3 chr1:1000000-2000000 --show 5`,
	Args: cobra.ExactArgs(2),
	RunE: func(cmd *cobra.Command, args []string) error {
		datasetPath := args[0]
		regionStr := args[1]

		// Parse region
		region, err := bams3.ParseRegion(regionStr)
		if err != nil {
			return fmt.Errorf("invalid region: %w", err)
		}

		// Open dataset
		reader, err := bams3.OpenDataset(datasetPath)
		if err != nil {
			return fmt.Errorf("failed to open dataset: %w", err)
		}
		defer reader.Close()

		// Query region
		fmt.Printf("Query: %s:%d-%d\n", region.Reference, region.Start, region.End)

		reads, err := reader.QueryRegion(region)
		if err != nil {
			return fmt.Errorf("query failed: %w", err)
		}

		fmt.Printf("Found %d reads in region\n", len(reads))

		if countOnly {
			return nil
		}

		// Show reads
		numToShow := showReads
		if numToShow > len(reads) {
			numToShow = len(reads)
		}

		if numToShow > 0 {
			fmt.Println()
			fmt.Printf("%-20s %12s %6s %s\n", "Read Name", "Position", "MapQ", "CIGAR")
			fmt.Println("------------------------------------------------------------")
			for i := 0; i < numToShow; i++ {
				read := reads[i]
				fmt.Printf("%-20s %12d %6d %s\n",
					read.Name,
					read.Position,
					read.MappingQuality,
					read.CIGAR)
			}
		}

		return nil
	},
}

func init() {
	queryCmd.Flags().BoolVar(&countOnly, "count", false,
		"Only show read count, don't display reads")
	queryCmd.Flags().IntVar(&showReads, "show", 10,
		"Number of reads to display (default 10, 0 for all)")
}
