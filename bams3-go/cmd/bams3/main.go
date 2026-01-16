package main

import (
	"fmt"
	"os"

	"github.com/spf13/cobra"
)

var rootCmd = &cobra.Command{
	Use:   "bams3",
	Short: "BAMS3 - Cloud-native alignment format tools",
	Long: `BAMS3 is a cloud-native alignment format designed for object storage.

This tool provides commands for converting BAM files to BAMS3 format,
querying BAMS3 datasets, and converting back to BAM format.`,
}

func main() {
	if err := rootCmd.Execute(); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
}

func init() {
	rootCmd.AddCommand(convertCmd)
	rootCmd.AddCommand(queryCmd)
	rootCmd.AddCommand(statsCmd)
	rootCmd.AddCommand(versionCmd)
}

var versionCmd = &cobra.Command{
	Use:   "version",
	Short: "Print version information",
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println("bams3-go version 0.1.0")
		fmt.Println("Cloud-native alignment format for S3")
	},
}
