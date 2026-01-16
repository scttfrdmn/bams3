#!/usr/bin/env python3
"""
Create a small synthetic BAM file for testing BAMS3 tools.
"""

try:
    import pysam
except ImportError:
    print("Error: pysam not installed")
    print("Install with: pip install pysam")
    exit(1)

import random
import sys

def create_test_bam(output_file, num_reads=1000):
    """Create a small test BAM file."""

    # Define header
    header = {
        'HD': {'VN': '1.6', 'SO': 'coordinate'},
        'SQ': [
            {'SN': 'chr1', 'LN': 10000000},
            {'SN': 'chr2', 'LN': 8000000},
            {'SN': 'chr3', 'LN': 6000000},
        ],
        'RG': [
            {
                'ID': 'rg1',
                'SM': 'test_sample',
                'LB': 'test_library',
                'PL': 'ILLUMINA'
            }
        ],
        'PG': [
            {
                'ID': 'test_generator',
                'VN': '1.0',
                'CL': 'create_test_bam.py'
            }
        ]
    }

    # Generate reads first, then sort
    bases = ['A', 'C', 'G', 'T']
    reads = []

    for i in range(num_reads):
        # Determine chromosome (distribution: 50% chr1, 30% chr2, 20% chr3)
        rand = random.random()
        if rand < 0.5:
            ref_id = 0  # chr1
            ref_length = 10000000
        elif rand < 0.8:
            ref_id = 1  # chr2
            ref_length = 8000000
        else:
            ref_id = 2  # chr3
            ref_length = 6000000

        # Generate read
        read = pysam.AlignedSegment()
        read.query_name = f'read_{i:06d}'
        read.flag = 0 if random.random() < 0.9 else 4  # 90% mapped, 10% unmapped

        if read.flag & 4:  # Unmapped
            read.reference_id = -1
            read.reference_start = -1
            read.mapping_quality = 0
            read.cigar = None
        else:  # Mapped
            read.reference_id = ref_id
            read.reference_start = random.randint(0, ref_length - 100)
            read.mapping_quality = random.randint(20, 60)
            read.cigar = [(0, 75)]  # 75M (match)

        # Generate sequence (75bp)
        sequence = ''.join(random.choice(bases) for _ in range(75))
        read.query_sequence = sequence

        # Generate quality scores (phred 20-40)
        qualities = [random.randint(20, 40) for _ in range(75)]
        read.query_qualities = qualities

        # Add tags
        read.set_tag('RG', 'rg1')

        reads.append(read)

    # Sort reads by reference ID and position
    reads.sort(key=lambda r: (r.reference_id if r.reference_id >= 0 else 999,
                              r.reference_start if r.reference_start >= 0 else 0))

    # Write sorted reads to BAM file
    bam = pysam.AlignmentFile(output_file, 'wb', header=header)
    for read in reads:
        bam.write(read)

    bam.close()

    print(f"✓ Created {output_file} with {num_reads} reads")

    # Create index
    print(f"Creating BAM index...")
    pysam.index(output_file)
    print(f"✓ Created {output_file}.bai")

if __name__ == '__main__':
    import os

    output_file = 'test_sample.bam'
    num_reads = 1000

    if len(sys.argv) > 1:
        output_file = sys.argv[1]
    if len(sys.argv) > 2:
        num_reads = int(sys.argv[2])

    print(f"Creating test BAM file: {output_file}")
    print(f"Number of reads: {num_reads}")
    print()

    create_test_bam(output_file, num_reads)

    # Show statistics
    print()
    print("File statistics:")
    bam = pysam.AlignmentFile(output_file, 'rb')

    ref_counts = {}
    total = 0
    mapped = 0

    for read in bam:
        total += 1
        if not read.is_unmapped:
            mapped += 1
            ref_name = bam.get_reference_name(read.reference_id)
            ref_counts[ref_name] = ref_counts.get(ref_name, 0) + 1

    print(f"  Total reads: {total}")
    print(f"  Mapped reads: {mapped}")
    print(f"  Unmapped reads: {total - mapped}")
    print()
    print("Reads per chromosome:")
    for ref, count in sorted(ref_counts.items()):
        print(f"  {ref}: {count}")

    bam.close()

    print()
    print("Test file ready!")
    print(f"  File: {output_file}")
    print(f"  Size: {os.path.getsize(output_file) / 1024:.1f} KB")
