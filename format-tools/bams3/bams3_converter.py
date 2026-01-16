#!/usr/bin/env python3
"""
BAMS3 Converter - Proof of Concept

Convert BAM files to BAMS3 (cloud-native alignment format).

Usage:
    python bams3_converter.py input.bam output.bams3 --chunk-size 1000000

This is a simplified proof-of-concept. Production version would:
- Use optimized binary formats
- Implement parallel processing
- Add comprehensive error handling
- Support streaming from S3
"""

import json
import sys
import argparse
from pathlib import Path
from datetime import datetime
from collections import defaultdict
from typing import Dict, List
import hashlib


def bam_to_bams3(
    input_bam: str,
    output_dir: str,
    chunk_size: int = 1000000,
    compression: str = 'none'
):
    """
    Convert BAM to BAMS3 format.

    Args:
        input_bam: Path to input BAM file
        output_dir: Output directory for BAMS3 dataset
        chunk_size: Genomic bases per chunk (default 1Mbp)
        compression: Compression algorithm (TODO)
    """
    try:
        import pysam
    except ImportError:
        print("Error: pysam is required. Install with: pip install pysam")
        sys.exit(1)

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    print(f"Converting {input_bam} to BAMS3 format...")
    print(f"Output directory: {output_dir}")
    print(f"Chunk size: {chunk_size:,} bp")

    # Open input BAM
    bam = pysam.AlignmentFile(input_bam, 'rb')

    # Create directory structure
    data_dir = output_path / 'data'
    index_dir = output_path / '_index'
    data_dir.mkdir(exist_ok=True)
    index_dir.mkdir(exist_ok=True)

    # Extract and save header
    print("\nExtracting header...")
    header = extract_header(bam)
    save_json(output_path / '_header.json', header)

    # Process reads into chunks
    print("\nProcessing reads into chunks...")
    chunks_metadata = []
    statistics = {
        'total_reads': 0,
        'mapped_reads': 0,
        'unmapped_reads': 0,
        'duplicate_reads': 0,
        'total_bases': 0
    }

    # Group reads by reference and position
    chunks = defaultdict(list)
    current_ref = None
    current_chunk_start = 0
    read_count = 0

    for read in bam:
        read_count += 1

        # Update statistics
        statistics['total_reads'] += 1
        if not read.is_unmapped:
            statistics['mapped_reads'] += 1
        else:
            statistics['unmapped_reads'] += 1
        if read.is_duplicate:
            statistics['duplicate_reads'] += 1
        if read.query_sequence:
            statistics['total_bases'] += len(read.query_sequence)

        # Determine chunk for this read
        if read.is_unmapped:
            chunk_key = ('unmapped', 0, 0)
        else:
            ref_name = bam.get_reference_name(read.reference_id)
            chunk_start = (read.reference_start // chunk_size) * chunk_size
            chunk_end = chunk_start + chunk_size
            chunk_key = (ref_name, chunk_start, chunk_end)

        chunks[chunk_key].append(read)

        if read_count % 100000 == 0:
            print(f"  Processed {read_count:,} reads...")

    print(f"\nTotal reads processed: {read_count:,}")
    print(f"Total chunks: {len(chunks):,}")

    # Write chunks to disk
    print("\nWriting chunks...")
    for idx, (chunk_key, reads) in enumerate(chunks.items()):
        ref_name, chunk_start, chunk_end = chunk_key

        # Create reference directory
        if ref_name == 'unmapped':
            chunk_dir = data_dir
            chunk_filename = 'unmapped.chunk'
        else:
            chunk_dir = data_dir / ref_name
            chunk_dir.mkdir(exist_ok=True)
            chunk_filename = f'{chunk_start:09d}-{chunk_end:09d}.chunk'

        chunk_path = chunk_dir / chunk_filename

        # Write chunk (simplified format for POC)
        chunk_data = write_chunk(chunk_path, reads)

        # Record chunk metadata
        chunks_metadata.append({
            'path': str(chunk_path.relative_to(output_path)),
            'reference': ref_name,
            'start': chunk_start,
            'end': chunk_end,
            'reads': len(reads),
            'size_bytes': chunk_path.stat().st_size,
            'compression': compression,
            'checksum': chunk_data['checksum'],
            'created': datetime.now().isoformat()
        })

        if (idx + 1) % 10 == 0:
            print(f"  Written {idx + 1}/{len(chunks)} chunks...")

    bam.close()

    # Calculate mean coverage (approximate)
    total_ref_length = sum(header['SQ'][i]['LN'] for i in range(len(header['SQ'])))
    if total_ref_length > 0:
        statistics['mean_coverage'] = statistics['total_bases'] / total_ref_length

    # Create metadata
    print("\nCreating metadata...")
    metadata = {
        'format': 'bams3',
        'version': '0.1.0-poc',
        'created': datetime.now().isoformat(),
        'created_by': 'bams3_converter.py POC',
        'source': {
            'file': input_bam,
            'format': 'BAM'
        },
        'statistics': statistics,
        'chunks': chunks_metadata,
        'compression': {
            'algorithm': compression,
        },
        'chunk_size': chunk_size,
    }

    save_json(output_path / '_metadata.json', metadata)

    # Create simple spatial index
    print("Creating spatial index...")
    create_spatial_index(output_path / '_index' / 'spatial.json', chunks_metadata)

    print("\n✓ Conversion complete!")
    print(f"\nDataset summary:")
    print(f"  Location: {output_path}")
    print(f"  Total reads: {statistics['total_reads']:,}")
    print(f"  Mapped reads: {statistics['mapped_reads']:,}")
    print(f"  Chunks: {len(chunks_metadata):,}")
    print(f"  Total size: {sum(c['size_bytes'] for c in chunks_metadata) / (1024**2):.1f} MB")

    print(f"\nTo query this dataset:")
    print(f"  python bams3_query.py {output_dir} chr1:1000000-2000000")


def extract_header(bam) -> Dict:
    """Extract SAM header as JSON."""
    header = bam.header.to_dict()

    # Simplify for readability
    return {
        'HD': header.get('HD', {}),
        'SQ': header.get('SQ', []),
        'RG': header.get('RG', []),
        'PG': header.get('PG', []),
    }


def write_chunk(chunk_path: Path, reads: List) -> Dict:
    """
    Write reads to chunk file.

    For POC, using simple JSON format. Production would use optimized binary.
    """
    chunk_data = []
    hasher = hashlib.sha256()

    for read in reads:
        # Simplified read record
        record = {
            'name': read.query_name,
            'flag': read.flag,
            'ref': read.reference_id,
            'pos': read.reference_start if not read.is_unmapped else -1,
            'mapq': read.mapping_quality,
            'cigar': read.cigarstring,
            'seq': read.query_sequence if read.query_sequence else '*',
            'qual': read.qual if read.qual else '*',
        }

        chunk_data.append(record)

        # Update checksum
        hasher.update(json.dumps(record, sort_keys=True).encode())

    # Write to file
    with open(chunk_path, 'w') as f:
        json.dump(chunk_data, f)

    return {
        'checksum': hasher.hexdigest()
    }


def create_spatial_index(index_path: Path, chunks_metadata: List[Dict]):
    """
    Create spatial index: position → chunk mapping.

    Simplified version for POC.
    """
    # Group chunks by reference
    by_reference = defaultdict(list)
    for chunk in chunks_metadata:
        ref = chunk['reference']
        by_reference[ref].append({
            'start': chunk['start'],
            'end': chunk['end'],
            'path': chunk['path'],
            'reads': chunk['reads']
        })

    # Sort each reference by position
    for ref in by_reference:
        by_reference[ref].sort(key=lambda x: x['start'])

    index = {
        'format': 'spatial_index',
        'version': '0.1.0',
        'references': dict(by_reference)
    }

    save_json(index_path, index)


def save_json(path: Path, data: Dict):
    """Save data as pretty-printed JSON."""
    with open(path, 'w') as f:
        json.dump(data, f, indent=2)


def main():
    parser = argparse.ArgumentParser(
        description='Convert BAM to BAMS3 format (Proof of Concept)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Convert local BAM
  python bams3_converter.py input.bam output.bams3

  # Custom chunk size
  python bams3_converter.py input.bam output.bams3 --chunk-size 5000000

Note: This is a proof-of-concept. Production implementation would:
  - Use optimized binary format (not JSON)
  - Support parallel processing
  - Handle large files efficiently
  - Support S3 input/output
  - Add compression (zstd, lz4)
        """
    )

    parser.add_argument('input_bam', help='Input BAM file')
    parser.add_argument('output_dir', help='Output BAMS3 directory')
    parser.add_argument(
        '--chunk-size',
        type=int,
        default=1000000,
        help='Chunk size in base pairs (default: 1000000)'
    )
    parser.add_argument(
        '--compression',
        default='none',
        choices=['none', 'zstd', 'lz4'],
        help='Compression algorithm (default: none, TODO)'
    )

    args = parser.parse_args()

    try:
        bam_to_bams3(
            args.input_bam,
            args.output_dir,
            args.chunk_size,
            args.compression
        )
    except Exception as e:
        print(f"\nError: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
