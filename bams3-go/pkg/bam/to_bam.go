package bam

import (
	"fmt"
	"os"
	"sort"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/scttfrdmn/bams3-go/pkg/bams3"
)

// ConvertBAMS3ToBAM converts a BAMS3 dataset back to BAM format
func ConvertBAMS3ToBAM(bams3Path string, outputBAM string) error {
	fmt.Printf("Converting BAMS3 to BAM format...\n")
	fmt.Printf("Input: %s\n", bams3Path)
	fmt.Printf("Output: %s\n\n", outputBAM)

	// Open BAMS3 dataset
	reader, err := bams3.OpenDataset(bams3Path)
	if err != nil {
		return fmt.Errorf("failed to open BAMS3 dataset: %w", err)
	}
	defer reader.Close()

	metadata := reader.GetMetadata()
	headerData := reader.GetHeader()

	fmt.Printf("Dataset: %s v%s\n", metadata.Format, metadata.Version)
	fmt.Printf("Total reads: %d\n", metadata.Statistics.TotalReads)
	fmt.Printf("Total chunks: %d\n\n", len(metadata.Chunks))

	// Convert header
	samHeader, err := convertToBamHeader(headerData)
	if err != nil {
		return fmt.Errorf("failed to convert header: %w", err)
	}

	// Create output BAM file
	outFile, err := os.Create(outputBAM)
	if err != nil {
		return fmt.Errorf("failed to create output file: %w", err)
	}
	defer outFile.Close()

	bamWriter, err := bam.NewWriter(outFile, samHeader, 1)
	if err != nil {
		return fmt.Errorf("failed to create BAM writer: %w", err)
	}
	defer bamWriter.Close()

	fmt.Println("Processing chunks...")

	// Sort chunks by reference and position
	chunks := make([]bams3.ChunkInfo, len(metadata.Chunks))
	copy(chunks, metadata.Chunks)
	sort.Slice(chunks, func(i, j int) bool {
		if chunks[i].Reference != chunks[j].Reference {
			return chunks[i].Reference < chunks[j].Reference
		}
		return chunks[i].Start < chunks[j].Start
	})

	totalReads := 0
	for i, chunk := range chunks {
		if (i+1)%5 == 0 || i == len(chunks)-1 {
			fmt.Printf("  Processing chunk %d/%d (%s:%d-%d)\r",
				i+1, len(chunks), chunk.Reference, chunk.Start, chunk.End)
		}

		// Load chunk from BAMS3
		reads, err := loadChunkReads(reader, chunk)
		if err != nil {
			return fmt.Errorf("failed to load chunk %s: %w", chunk.Path, err)
		}

		// Convert and write each read
		for _, read := range reads {
			samRecord, err := convertToBamRecord(read, samHeader)
			if err != nil {
				return fmt.Errorf("failed to convert read %s: %w", read.Name, err)
			}

			if err := bamWriter.Write(samRecord); err != nil {
				return fmt.Errorf("failed to write read: %w", err)
			}

			totalReads++
		}
	}

	fmt.Println()
	fmt.Printf("\n✓ Conversion complete!\n")
	fmt.Printf("  Reads written: %d\n", totalReads)
	fmt.Printf("  Output: %s\n\n", outputBAM)

	// Close writer to flush
	if err := bamWriter.Close(); err != nil {
		return fmt.Errorf("failed to close BAM writer: %w", err)
	}

	// Create index
	fmt.Println("Creating BAM index...")
	if err := indexBAM(outputBAM); err != nil {
		fmt.Printf("Warning: Could not create index: %v\n", err)
		fmt.Println("  (This is normal if reads are not coordinate-sorted)")
	} else {
		fmt.Printf("✓ Index created: %s.bai\n", outputBAM)
	}

	return nil
}

// loadChunkReads loads all reads from a specific chunk
func loadChunkReads(reader *bams3.Reader, chunk bams3.ChunkInfo) ([]bams3.Read, error) {
	// Query the entire chunk region
	region := bams3.Region{
		Reference: chunk.Reference,
		Start:     chunk.Start,
		End:       chunk.End,
	}

	return reader.QueryRegion(region)
}

// convertToBamHeader converts BAMS3 header to sam.Header
func convertToBamHeader(headerData bams3.Header) (*sam.Header, error) {
	// Create references first
	var refs []*sam.Reference
	for _, sq := range headerData.SQ {
		var length int
		if ln, ok := sq["LN"]; ok {
			fmt.Sscanf(ln, "%d", &length)
		}

		name := sq["SN"]
		ref, err := sam.NewReference(name, "", "", length, nil, nil)
		if err != nil {
			return nil, fmt.Errorf("failed to create reference %s: %w", name, err)
		}
		refs = append(refs, ref)
	}

	// Create header with references
	header, err := sam.NewHeader(nil, refs)
	if err != nil {
		return nil, fmt.Errorf("failed to create header: %w", err)
	}

	// Set version
	if vn, ok := headerData.HD["VN"]; ok {
		header.Version = vn
	}

	// Set sort order
	if so, ok := headerData.HD["SO"]; ok {
		switch so {
		case "coordinate":
			header.SortOrder = sam.Coordinate
		case "queryname":
			header.SortOrder = sam.QueryName
		case "unsorted":
			header.SortOrder = sam.Unsorted
		default:
			header.SortOrder = sam.UnknownOrder
		}
	}

	// Add read groups (simplified)
	for _, rg := range headerData.RG {
		if id, ok := rg["ID"]; ok {
			// Create minimal read group
			readGroup := &sam.ReadGroup{}
			readGroup.SetName(id)
			header.AddReadGroup(readGroup)
		}
	}

	// Add comments
	header.Comments = headerData.CO

	return header, nil
}

// convertToBamRecord converts BAMS3 read to sam.Record
func convertToBamRecord(read bams3.Read, header *sam.Header) (*sam.Record, error) {
	record := &sam.Record{
		Name:  read.Name,
		Flags: sam.Flags(read.Flag),
		MapQ:  read.MappingQuality,
	}

	// Set reference
	if read.ReferenceID >= 0 && read.ReferenceID < len(header.Refs()) {
		record.Ref = header.Refs()[read.ReferenceID]
		record.Pos = read.Position
	} else {
		record.Ref = nil
		record.Pos = -1
	}

	// Parse CIGAR
	if read.CIGAR != "" && read.CIGAR != "*" {
		cigar, err := sam.ParseCigar([]byte(read.CIGAR))
		if err != nil {
			return nil, fmt.Errorf("failed to parse CIGAR: %w", err)
		}
		record.Cigar = cigar
	}

	// Set sequence
	if read.Sequence != "" && read.Sequence != "*" {
		seq := sam.NewSeq([]byte(read.Sequence))
		record.Seq = seq
	}

	// Set quality
	if read.Quality != "" && read.Quality != "*" {
		qual := make([]byte, len(read.Quality))
		for i, q := range read.Quality {
			qual[i] = byte(q) - 33
		}
		record.Qual = qual
	}

	// Mate information (unmapped for POC)
	record.MateRef = nil
	record.MatePos = -1
	record.TempLen = 0

	// Tags (if present)
	if len(read.Tags) > 0 {
		for tagName, tagValue := range read.Tags {
			if len(tagName) == 2 {
				// Create auxiliary field
				// This is simplified - full implementation needs proper tag type handling
				aux := sam.Aux([]byte(fmt.Sprintf("%s:Z:%v", tagName, tagValue)))
				record.AuxFields = append(record.AuxFields, aux)
			}
		}
	}

	return record, nil
}

// indexBAM creates a BAM index file using samtools
func indexBAM(bamPath string) error {
	// For now, we'll skip automatic indexing
	// Users can run: samtools index output.bam
	// Future: implement proper indexing or call samtools
	return fmt.Errorf("automatic indexing not yet implemented - run 'samtools index %s'", bamPath)
}
