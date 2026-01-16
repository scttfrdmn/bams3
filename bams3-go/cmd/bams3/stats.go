package main

import (
	"fmt"

	"github.com/scttfrdmn/bams3-go/pkg/bams3"
	"github.com/spf13/cobra"
)

var statsCmd = &cobra.Command{
	Use:   "stats <dataset.bams3>",
	Short: "Show statistics for a BAMS3 dataset",
	Long: `Display statistics for a BAMS3 dataset.

Statistics are read instantly from the metadata file without
scanning the entire dataset.

Example:
  bams3 stats sample.bams3`,
	Args: cobra.ExactArgs(1),
	RunE: func(cmd *cobra.Command, args []string) error {
		datasetPath := args[0]

		// Open dataset
		reader, err := bams3.OpenDataset(datasetPath)
		if err != nil {
			return fmt.Errorf("failed to open dataset: %w", err)
		}
		defer reader.Close()

		// Get metadata and statistics
		metadata := reader.GetMetadata()
		stats := reader.GetStatistics()
		header := reader.GetHeader()

		// Display information
		fmt.Println("===========================================")
		fmt.Println("BAMS3 Dataset Statistics")
		fmt.Println("===========================================")
		fmt.Println()
		fmt.Printf("Format: %s v%s\n", metadata.Format, metadata.Version)
		fmt.Printf("Created: %s\n", metadata.Created.Format("2006-01-02 15:04:05"))
		fmt.Printf("Created by: %s\n", metadata.CreatedBy)
		if metadata.Source.File != "" {
			fmt.Printf("Source: %s (%s)\n", metadata.Source.File, metadata.Source.Format)
		}
		fmt.Println()

		fmt.Println("Statistics:")
		fmt.Printf("  Total reads: %d\n", stats.TotalReads)
		fmt.Printf("  Mapped reads: %d (%.2f%%)\n", stats.MappedReads,
			float64(stats.MappedReads)/float64(stats.TotalReads)*100)
		fmt.Printf("  Unmapped reads: %d (%.2f%%)\n", stats.UnmappedReads,
			float64(stats.UnmappedReads)/float64(stats.TotalReads)*100)
		fmt.Printf("  Duplicate reads: %d\n", stats.DuplicateReads)
		fmt.Printf("  Total bases: %d\n", stats.TotalBases)
		fmt.Printf("  Mean coverage: %.2fx\n", stats.MeanCoverage)
		fmt.Println()

		fmt.Println("Structure:")
		fmt.Printf("  Chunk size: %d bp\n", metadata.ChunkSize)
		fmt.Printf("  Total chunks: %d\n", len(metadata.Chunks))
		fmt.Printf("  Compression: %s\n", metadata.Compression.Algorithm)
		fmt.Println()

		fmt.Println("References:")
		for _, sq := range header.SQ {
			fmt.Printf("  %s: %s bp\n", sq["SN"], sq["LN"])
		}

		return nil
	},
}
