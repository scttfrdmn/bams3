#!/usr/bin/env python3
"""
BAMS3 Query Tool - Proof of Concept

Query BAMS3 datasets by genomic region.

Usage:
    python bams3_query.py sample.bams3 chr1:1000000-2000000
    python bams3_query.py sample.bams3 chr1  # Entire chromosome
    python bams3_query.py sample.bams3 --stats  # Show statistics

Demonstrates efficient querying - only downloads relevant chunks.
"""

import json
import sys
import argparse
from pathlib import Path
from typing import Dict, List, Tuple


class BAMS3Dataset:
    """Simple BAMS3 dataset reader."""

    def __init__(self, dataset_path: str):
        self.path = Path(dataset_path)

        if not self.path.exists():
            raise FileNotFoundError(f"Dataset not found: {dataset_path}")

        # Load metadata
        metadata_file = self.path / '_metadata.json'
        with open(metadata_file) as f:
            self.metadata = json.load(f)

        # Load header
        header_file = self.path / '_header.json'
        with open(header_file) as f:
            self.header = json.load(f)

        # Load spatial index
        index_file = self.path / '_index' / 'spatial.json'
        with open(index_file) as f:
            self.index = json.load(f)

    def query_region(self, reference: str, start: int, end: int):
        """
        Query reads in genomic region.

        Only loads chunks that overlap the query region.
        """
        # Find overlapping chunks using index
        chunks = self._find_overlapping_chunks(reference, start, end)

        print(f"\nQuery: {reference}:{start:,}-{end:,}")
        print(f"Chunks to load: {len(chunks)}")

        total_reads = 0
        matching_reads = 0

        for chunk_info in chunks:
            chunk_path = self.path / chunk_info['path']
            print(f"  Loading chunk: {chunk_info['path']} ({chunk_info['reads']} reads)")

            # Load chunk
            with open(chunk_path) as f:
                reads = json.load(f)

            total_reads += len(reads)

            # Filter to exact region
            for read in reads:
                if read['pos'] >= start and read['pos'] < end:
                    matching_reads += 1
                    yield read

        print(f"\nScanned {total_reads:,} reads from chunks")
        print(f"Found {matching_reads:,} reads in exact region")

    def query_chromosome(self, reference: str):
        """Query all reads on a chromosome."""
        if reference not in self.index['references']:
            print(f"Error: Reference '{reference}' not found")
            return

        chunks = self.index['references'][reference]
        print(f"\nQuerying entire chromosome: {reference}")
        print(f"Chunks: {len(chunks)}")

        total_reads = 0

        for chunk_info in chunks:
            chunk_path = self.path / chunk_info['path']
            print(f"  Loading chunk: {chunk_info['path']}")

            with open(chunk_path) as f:
                reads = json.load(f)

            total_reads += len(reads)

            for read in reads:
                yield read

        print(f"\nTotal reads: {total_reads:,}")

    def _find_overlapping_chunks(self, reference: str, start: int, end: int) -> List[Dict]:
        """Find chunks that overlap query region."""
        if reference not in self.index['references']:
            return []

        overlapping = []

        for chunk in self.index['references'][reference]:
            # Check if chunk overlaps query region
            chunk_start = chunk['start']
            chunk_end = chunk['end']

            if chunk_end > start and chunk_start < end:
                overlapping.append(chunk)

        return overlapping

    def get_statistics(self) -> Dict:
        """Get dataset statistics."""
        return self.metadata['statistics']

    def list_references(self) -> List[str]:
        """List available references."""
        return list(self.index['references'].keys())

    def get_info(self) -> Dict:
        """Get dataset information."""
        return {
            'format': self.metadata['format'],
            'version': self.metadata['version'],
            'created': self.metadata['created'],
            'chunks': len(self.metadata['chunks']),
            'chunk_size': self.metadata['chunk_size'],
            'references': self.list_references(),
            'statistics': self.get_statistics()
        }


def parse_region(region_str: str) -> Tuple[str, int, int]:
    """
    Parse region string like 'chr1:1000000-2000000'.

    Returns:
        (reference, start, end)
    """
    if ':' not in region_str:
        # Just chromosome name
        return region_str, None, None

    ref, pos = region_str.split(':')

    if '-' not in pos:
        raise ValueError(f"Invalid region format: {region_str}")

    start, end = pos.split('-')
    start = int(start.replace(',', ''))
    end = int(end.replace(',', ''))

    return ref, start, end


def print_read(read: Dict, format: str = 'simple'):
    """Print read in various formats."""
    if format == 'simple':
        print(f"{read['name']}\t{read['pos']:,}\t{read['mapq']}\t{read['cigar']}")
    elif format == 'sam':
        print(f"{read['name']}\t{read['flag']}\t{read['ref']}\t{read['pos']}\t"
              f"{read['mapq']}\t{read['cigar']}\t*\t0\t0\t{read['seq']}\t{read['qual']}")


def main():
    parser = argparse.ArgumentParser(
        description='Query BAMS3 datasets',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Query specific region
  python bams3_query.py sample.bams3 chr1:1000000-2000000

  # Query entire chromosome
  python bams3_query.py sample.bams3 chr1

  # Show dataset statistics
  python bams3_query.py sample.bams3 --stats

  # Show dataset info
  python bams3_query.py sample.bams3 --info

  # Count reads in region
  python bams3_query.py sample.bams3 chr1:1000000-2000000 --count

  # Output in SAM format
  python bams3_query.py sample.bams3 chr1:1000000-2000000 --format sam
        """
    )

    parser.add_argument('dataset', help='BAMS3 dataset directory')
    parser.add_argument('region', nargs='?', help='Region to query (chr:start-end)')

    parser.add_argument('--stats', action='store_true', help='Show statistics')
    parser.add_argument('--info', action='store_true', help='Show dataset info')
    parser.add_argument('--count', action='store_true', help='Count reads only')
    parser.add_argument('--format', default='simple', choices=['simple', 'sam'],
                        help='Output format')
    parser.add_argument('--head', type=int, help='Show only first N reads')

    args = parser.parse_args()

    try:
        # Open dataset
        dataset = BAMS3Dataset(args.dataset)

        # Show statistics
        if args.stats:
            stats = dataset.get_statistics()
            print("\nDataset Statistics:")
            print("=" * 50)
            for key, value in stats.items():
                if isinstance(value, float):
                    print(f"  {key}: {value:.2f}")
                else:
                    print(f"  {key}: {value:,}")
            return

        # Show info
        if args.info:
            info = dataset.get_info()
            print("\nDataset Information:")
            print("=" * 50)
            print(f"  Format: {info['format']} v{info['version']}")
            print(f"  Created: {info['created']}")
            print(f"  Chunks: {info['chunks']:,}")
            print(f"  Chunk size: {info['chunk_size']:,} bp")
            print(f"  References: {', '.join(info['references'][:10])}")
            if len(info['references']) > 10:
                print(f"    ... and {len(info['references']) - 10} more")
            print()
            print("Statistics:")
            for key, value in info['statistics'].items():
                if isinstance(value, float):
                    print(f"  {key}: {value:.2f}")
                else:
                    print(f"  {key}: {value:,}")
            return

        # Query region
        if not args.region:
            print("Error: region required (or use --stats, --info)")
            sys.exit(1)

        # Parse region
        ref, start, end = parse_region(args.region)

        # Query
        if start is None:
            # Entire chromosome
            reads = dataset.query_chromosome(ref)
        else:
            # Specific region
            reads = dataset.query_region(ref, start, end)

        # Output reads
        count = 0
        for read in reads:
            count += 1

            if args.count:
                continue

            if args.head and count > args.head:
                print(f"\n(Showing first {args.head} reads)")
                break

            print_read(read, args.format)

        if args.count:
            print(f"\nTotal reads: {count:,}")

    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
