#!/usr/bin/env python3
"""
Convert BAMS3 format back to standard BAM format.
This enables bidirectional conversion for compatibility.
"""

import json
import sys
import os
from pathlib import Path
import pysam


def bams3_to_bam(bams3_dir, output_bam):
    """
    Convert BAMS3 dataset back to standard BAM format.

    Args:
        bams3_dir: Path to BAMS3 directory
        output_bam: Output BAM file path
    """
    bams3_path = Path(bams3_dir)

    # Read metadata
    metadata_file = bams3_path / '_metadata.json'
    if not metadata_file.exists():
        raise ValueError(f"Metadata file not found: {metadata_file}")

    with open(metadata_file) as f:
        metadata = json.load(f)

    print(f"Converting BAMS3 → BAM: {bams3_dir} → {output_bam}")
    print(f"Format: {metadata['format']} v{metadata['version']}")
    print(f"Total reads: {metadata['statistics']['total_reads']:,}")
    print(f"Total chunks: {len(metadata['chunks'])}")
    print()

    # Read header
    header_file = bams3_path / '_header.json'
    if not header_file.exists():
        raise ValueError(f"Header file not found: {header_file}")

    with open(header_file) as f:
        header_data = json.load(f)

    # Reconstruct header for pysam
    header = {
        'HD': header_data['HD'],
        'SQ': header_data['SQ']
    }
    if 'PG' in header_data:
        header['PG'] = header_data['PG']
    if 'RG' in header_data:
        header['RG'] = header_data['RG']
    if 'CO' in header_data:
        header['CO'] = header_data['CO']

    # Open output BAM
    outfile = pysam.AlignmentFile(output_bam, 'wb', header=header)

    reads_written = 0

    # Process chunks in order (by reference and position)
    sorted_chunks = sorted(metadata['chunks'],
                          key=lambda c: (c.get('reference', 'zzz'),
                                       c.get('start', 0)))

    for i, chunk_info in enumerate(sorted_chunks):
        chunk_path = bams3_path / chunk_info['path']

        if not chunk_path.exists():
            print(f"Warning: Chunk not found: {chunk_path}", file=sys.stderr)
            continue

        # Load chunk
        with open(chunk_path) as f:
            chunk_data = json.load(f)

        # Convert each read back to pysam format
        for read_dict in chunk_data:
            read = pysam.AlignedSegment()

            # Required fields
            read.query_name = read_dict['name']
            read.flag = read_dict['flag']

            # Reference and position
            ref_id = read_dict.get('ref', -1)
            if ref_id >= 0:
                read.reference_id = ref_id
                read.reference_start = read_dict['pos']
            else:
                read.reference_id = -1
                read.reference_start = -1

            read.mapping_quality = read_dict['mapq']

            # CIGAR string
            cigar = read_dict.get('cigar')
            if cigar and cigar != '*':
                read.cigarstring = cigar

            # Mate info (set to unmapped for POC format which doesn't store mate info)
            read.next_reference_id = -1
            read.next_reference_start = -1
            read.template_length = 0

            # Sequence and quality
            seq = read_dict.get('seq', '*')
            if seq and seq != '*':
                read.query_sequence = seq

            qual = read_dict.get('qual', '*')
            if qual and qual != '*':
                read.query_qualities = pysam.qualitystring_to_array(qual)

            # Tags (if present in future versions)
            if 'tags' in read_dict:
                for tag_name, tag_value in read_dict['tags'].items():
                    # Infer tag type from value
                    if isinstance(tag_value, int):
                        read.set_tag(tag_name, tag_value, value_type='i')
                    elif isinstance(tag_value, float):
                        read.set_tag(tag_name, tag_value, value_type='f')
                    elif isinstance(tag_value, str):
                        read.set_tag(tag_name, tag_value, value_type='Z')

            outfile.write(read)
            reads_written += 1

        if (i + 1) % 5 == 0 or i == len(sorted_chunks) - 1:
            print(f"Processed {i + 1}/{len(sorted_chunks)} chunks ({reads_written:,} reads)", end='\r')

    print()  # New line after progress

    outfile.close()

    print(f"\n✓ Conversion complete!")
    print(f"  Reads written: {reads_written:,}")
    print(f"  Output: {output_bam}")

    # Create index
    print("\nCreating BAM index...")
    try:
        pysam.index(output_bam)
        print(f"✓ Index created: {output_bam}.bai")
    except Exception as e:
        print(f"Warning: Could not create index: {e}", file=sys.stderr)
        print("  (Reads may not be properly sorted)", file=sys.stderr)

    return reads_written


def main():
    if len(sys.argv) != 3:
        print("Usage: bams3_to_bam.py <input.bams3/> <output.bam>")
        print()
        print("Convert BAMS3 format back to standard BAM format.")
        print()
        print("Example:")
        print("  python3 bams3_to_bam.py test_sample.bams3 reconstructed.bam")
        sys.exit(1)

    bams3_dir = sys.argv[1]
    output_bam = sys.argv[2]

    if not os.path.isdir(bams3_dir):
        print(f"Error: BAMS3 directory not found: {bams3_dir}", file=sys.stderr)
        sys.exit(1)

    try:
        bams3_to_bam(bams3_dir, output_bam)
    except Exception as e:
        print(f"\nError: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
