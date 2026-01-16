package bam

import (
	"fmt"
	"os"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/scttfrdmn/bams3-go/pkg/bams3"
)

// ConvertBAMToBAMS3 converts a BAM file to BAMS3 format
func ConvertBAMToBAMS3(bamPath string, outputPath string, chunkSize int, compression string, format string) error {
	fmt.Printf("Converting %s to BAMS3 format...\n", bamPath)
	fmt.Printf("Output directory: %s\n", outputPath)
	fmt.Printf("Chunk size: %d bp\n", chunkSize)
	fmt.Printf("Format: %s\n", format)
	fmt.Printf("Compression: %s\n\n", compression)

	// Open BAM file
	f, err := os.Open(bamPath)
	if err != nil {
		return fmt.Errorf("failed to open BAM file: %w", err)
	}
	defer f.Close()

	bamFile, err := bam.NewReader(f, 1)
	if err != nil {
		return fmt.Errorf("failed to create BAM reader: %w", err)
	}
	defer bamFile.Close()

	// Create BAMS3 writer
	writer, err := bams3.NewWriter(outputPath, chunkSize, compression, format)
	if err != nil {
		return fmt.Errorf("failed to create BAMS3 writer: %w", err)
	}

	// Convert and set header
	header := convertHeader(bamFile.Header())
	writer.SetHeader(header)

	// Set source
	writer.SetSource(bams3.Source{
		File:   bamPath,
		Format: "BAM",
	})

	fmt.Println("Processing reads into chunks...")

	// Process reads
	readCount := 0
	for {
		record, err := bamFile.Read()
		if err != nil {
			if err.Error() == "EOF" {
				break
			}
			return fmt.Errorf("failed to read BAM record: %w", err)
		}

		readCount++
		if readCount%100000 == 0 {
			fmt.Printf("  Processed %d reads...\n", readCount)
		}

		// Convert BAM record to BAMS3 Read
		read := convertRecord(record)

		// Get reference name
		refName := "*"
		if record.Ref != nil {
			refName = record.Ref.Name()
		}

		// Add to writer
		if err := writer.AddRead(read, refName); err != nil {
			return fmt.Errorf("failed to add read: %w", err)
		}
	}

	fmt.Printf("\nTotal reads processed: %d\n\n", readCount)

	// Finalize (writes all chunks and metadata)
	return writer.Finalize()
}

// convertHeader converts sam.Header to bams3.Header
func convertHeader(samHeader *sam.Header) bams3.Header {
	header := bams3.Header{
		HD: make(map[string]string),
		SQ: []map[string]string{},
		RG: []map[string]string{},
		PG: []map[string]string{},
		CO: []string{},
	}

	// Convert HD
	if samHeader.Version != "" {
		header.HD["VN"] = samHeader.Version
	}
	if samHeader.SortOrder != sam.UnknownOrder {
		header.HD["SO"] = samHeader.SortOrder.String()
	}

	// Convert SQ (reference sequences)
	for _, ref := range samHeader.Refs() {
		sq := map[string]string{
			"SN": ref.Name(),
			"LN": fmt.Sprintf("%d", ref.Len()),
		}
		header.SQ = append(header.SQ, sq)
	}

	// Convert RG (read groups)
	for _, rg := range samHeader.RGs() {
		rgMap := map[string]string{
			"ID": rg.Name(),
		}
		// Add other RG fields if available
		header.RG = append(header.RG, rgMap)
	}

	// Convert PG (programs)
	for _, pg := range samHeader.Progs() {
		pgMap := map[string]string{
			"ID": fmt.Sprintf("%d", pg.ID()),
		}
		if pg.Name() != "" {
			pgMap["PN"] = pg.Name()
		}
		if pg.Version() != "" {
			pgMap["VN"] = pg.Version()
		}
		// CommandLine not available in this version of biogo/hts
		header.PG = append(header.PG, pgMap)
	}

	// Convert CO (comments)
	header.CO = samHeader.Comments

	return header
}

// convertRecord converts a sam.Record to bams3.Read
func convertRecord(record *sam.Record) bams3.Read {
	read := bams3.Read{
		Name:           record.Name,
		Flag:           int(record.Flags),
		Position:       record.Pos,
		MappingQuality: uint8(record.MapQ),
		CIGAR:          record.Cigar.String(),
		Sequence:       string(record.Seq.Expand()),
		Tags:           make(map[string]interface{}),
	}

	// Convert CIGAR to operations (for binary format)
	if record.Cigar != nil {
		read.CIGAROps = make([]bams3.CIGAROperation, len(record.Cigar))
		for i, op := range record.Cigar {
			read.CIGAROps[i] = bams3.CIGAROperation{
				Type:   byte(op.Type()),
				Length: op.Len(),
			}
		}
	}

	// Set reference ID
	if record.Ref != nil {
		read.ReferenceID = record.Ref.ID()
	} else {
		read.ReferenceID = -1
	}

	// Mate information
	if record.MateRef != nil {
		read.MateReferenceID = record.MateRef.ID()
	} else {
		read.MateReferenceID = -1
	}
	read.MatePosition = record.MatePos
	read.TemplateLength = record.TempLen

	// Quality scores
	if record.Qual != nil {
		qual := make([]byte, len(record.Qual))
		for i, q := range record.Qual {
			qual[i] = byte(q + 33) // Convert to ASCII
		}
		read.Quality = string(qual)
	} else {
		read.Quality = "*"
	}

	// Convert tags
	for _, tag := range record.AuxFields {
		tagStr := tag.String()
		// Simple tag parsing - can be improved
		if len(tagStr) >= 2 {
			read.Tags[tagStr[:2]] = tagStr[5:]
		}
	}

	return read
}
