package bam

import (
	"fmt"
	"io"
	"os"
	"sort"
	"strings"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/scttfrdmn/bams3-go/pkg/bams3"
)

// ToBamOptions configures BAM export behavior
type ToBamOptions struct {
	Region      string // Optional region filter (chr:start-end or chr)
	CreateIndex bool   // Whether to create BAM index
}

// ConvertBAMS3ToBAM converts a BAMS3 dataset back to BAM format
func ConvertBAMS3ToBAM(bams3Path string, outputBAM string, opts ToBamOptions) error {
	// Determine if streaming to stdout
	isStdout := (outputBAM == "-")

	if !isStdout {
		fmt.Printf("Converting BAMS3 to BAM format...\n")
		fmt.Printf("Input: %s\n", bams3Path)
		if opts.Region != "" {
			fmt.Printf("Region: %s\n", opts.Region)
		}
		fmt.Printf("Output: %s\n\n", outputBAM)
	}

	// Open BAMS3 dataset
	reader, err := bams3.OpenDataset(bams3Path)
	if err != nil {
		return fmt.Errorf("failed to open BAMS3 dataset: %w", err)
	}
	defer reader.Close()

	metadata := reader.GetMetadata()
	headerData := reader.GetHeader()

	if !isStdout {
		fmt.Printf("Dataset: %s v%s\n", metadata.Format, metadata.Version)
		fmt.Printf("Total reads: %d\n", metadata.Statistics.TotalReads)
		fmt.Printf("Total chunks: %d\n\n", len(metadata.Chunks))
	}

	// Convert header
	samHeader, err := convertToBamHeader(headerData)
	if err != nil {
		return fmt.Errorf("failed to convert header: %w", err)
	}

	// Create output stream
	var outWriter io.Writer
	var outFile *os.File

	if isStdout {
		outWriter = os.Stdout
	} else {
		outFile, err = os.Create(outputBAM)
		if err != nil {
			return fmt.Errorf("failed to create output file: %w", err)
		}
		defer outFile.Close()
		outWriter = outFile
	}

	bamWriter, err := bam.NewWriter(outWriter, samHeader, 1)
	if err != nil {
		return fmt.Errorf("failed to create BAM writer: %w", err)
	}
	defer bamWriter.Close()

	// Determine which chunks to process
	var chunksToProcess []bams3.ChunkInfo

	if opts.Region != "" {
		// Parse region and filter chunks
		region, err := parseRegion(opts.Region)
		if err != nil {
			return fmt.Errorf("invalid region: %w", err)
		}

		// Filter chunks that overlap the region
		for _, chunk := range metadata.Chunks {
			if overlapsRegion(chunk, region) {
				chunksToProcess = append(chunksToProcess, chunk)
			}
		}

		if !isStdout {
			fmt.Printf("Processing %d chunks overlapping region...\n", len(chunksToProcess))
		}
	} else {
		// Process all chunks
		chunksToProcess = metadata.Chunks
		if !isStdout {
			fmt.Println("Processing chunks...")
		}
	}

	// Sort chunks by reference and position
	sort.Slice(chunksToProcess, func(i, j int) bool {
		if chunksToProcess[i].Reference != chunksToProcess[j].Reference {
			return chunksToProcess[i].Reference < chunksToProcess[j].Reference
		}
		return chunksToProcess[i].Start < chunksToProcess[j].Start
	})

	totalReads := 0
	for i, chunk := range chunksToProcess {
		if !isStdout && ((i+1)%5 == 0 || i == len(chunksToProcess)-1) {
			fmt.Printf("  Processing chunk %d/%d (%s:%d-%d)\r",
				i+1, len(chunksToProcess), chunk.Reference, chunk.Start, chunk.End)
		}

		// Load chunk from BAMS3
		var reads []bams3.Read
		if opts.Region != "" {
			// Use region query to filter reads
			region, _ := parseRegion(opts.Region)
			reads, err = reader.QueryRegion(region)
			if err != nil {
				return fmt.Errorf("failed to query region: %w", err)
			}
		} else {
			// Load entire chunk
			reads, err = loadChunkReads(reader, chunk)
			if err != nil {
				return fmt.Errorf("failed to load chunk %s: %w", chunk.Path, err)
			}
		}

		// Convert and write each read
		for _, read := range reads {
			// Apply region filter if specified
			if opts.Region != "" {
				region, _ := parseRegion(opts.Region)
				if !readOverlapsRegion(read, region) {
					continue
				}
			}

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

	// Close writer to flush
	if err := bamWriter.Close(); err != nil {
		return fmt.Errorf("failed to close BAM writer: %w", err)
	}

	if !isStdout {
		fmt.Println()
		fmt.Printf("\n✓ Conversion complete!\n")
		fmt.Printf("  Reads written: %d\n", totalReads)
		fmt.Printf("  Output: %s\n\n", outputBAM)

		// Create index if requested and not stdout
		if opts.CreateIndex {
			fmt.Println("Creating BAM index...")
			if err := indexBAM(outputBAM); err != nil {
				fmt.Printf("Warning: Could not create index: %v\n", err)
				fmt.Println("  Run 'samtools index %s' to create index manually\n", outputBAM)
			} else {
				fmt.Printf("✓ Index created: %s.bai\n", outputBAM)
			}
		}
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

// parseRegion parses a region string (chr:start-end or chr) into a Region struct
func parseRegion(regionStr string) (bams3.Region, error) {
	parts := strings.Split(regionStr, ":")
	if len(parts) == 1 {
		// Just chromosome name
		return bams3.Region{
			Reference: parts[0],
			Start:     0,
			End:       -1, // Entire chromosome
		}, nil
	}

	if len(parts) != 2 {
		return bams3.Region{}, fmt.Errorf("invalid region format: %s (expected chr:start-end or chr)", regionStr)
	}

	ref := parts[0]
	coords := strings.Split(parts[1], "-")
	if len(coords) != 2 {
		return bams3.Region{}, fmt.Errorf("invalid coordinates: %s (expected start-end)", parts[1])
	}

	var start, end int
	if _, err := fmt.Sscanf(coords[0], "%d", &start); err != nil {
		return bams3.Region{}, fmt.Errorf("invalid start coordinate: %s", coords[0])
	}
	if _, err := fmt.Sscanf(coords[1], "%d", &end); err != nil {
		return bams3.Region{}, fmt.Errorf("invalid end coordinate: %s", coords[1])
	}

	return bams3.Region{
		Reference: ref,
		Start:     start,
		End:       end,
	}, nil
}

// overlapsRegion checks if a chunk overlaps with a query region
func overlapsRegion(chunk bams3.ChunkInfo, region bams3.Region) bool {
	if chunk.Reference != region.Reference {
		return false
	}

	// If region is entire chromosome (end = -1), include chunk
	if region.End == -1 {
		return true
	}

	// Check for overlap: chunk.Start < region.End && chunk.End > region.Start
	return chunk.Start < region.End && chunk.End > region.Start
}

// readOverlapsRegion checks if a read overlaps with a query region
func readOverlapsRegion(read bams3.Read, region bams3.Region) bool {
	// Get reference name from read's reference ID
	// This is a simplified check - assumes reference names match
	// A full implementation would map reference IDs to names

	// If region is entire chromosome, include all reads from that chromosome
	if region.End == -1 {
		return true
	}

	// Calculate read end position from CIGAR
	readEnd := read.Position + calculateReadLength(read.CIGAR)

	// Check for overlap
	return read.Position < region.End && readEnd > region.Start
}

// calculateReadLength calculates the reference span of a read from its CIGAR
func calculateReadLength(cigar string) int {
	if cigar == "" || cigar == "*" {
		return 0
	}

	length := 0
	currentNum := 0

	for _, ch := range cigar {
		if ch >= '0' && ch <= '9' {
			currentNum = currentNum*10 + int(ch-'0')
		} else {
			// CIGAR operation
			switch ch {
			case 'M', 'D', 'N', '=', 'X':
				// These operations consume reference
				length += currentNum
			}
			currentNum = 0
		}
	}

	return length
}
