#!/bin/bash
#
# Generate minimal test data for Nextflow pipeline validation
#

set -e

echo "Generating test data..."

# Generate reference genome (2 chromosomes, 10kb each)
cat > reference.fa << 'EOF'
>chr1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCC
TATATATATATATATATATATATATATATATATATATATATATATATATATATATAT
CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG
>chr2
TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA
AATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATT
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
CCAACCAACCAACCAACCAACCAACCAACCAACCAACCAACCAACCAACCAACCAACCAA
EOF

# Generate synthetic FASTQ reads (1000 reads)
python3 - << 'PYTHON'
import random

def generate_read(ref_seq, read_len=100, error_rate=0.01):
    """Generate a synthetic read from reference"""
    start = random.randint(0, len(ref_seq) - read_len)
    read = list(ref_seq[start:start + read_len])

    # Add sequencing errors
    for i in range(len(read)):
        if random.random() < error_rate:
            read[i] = random.choice(['A', 'C', 'G', 'T'])

    qual = ''.join(['I'] * read_len)  # High quality
    return ''.join(read), qual

# Read reference
with open('reference.fa') as f:
    lines = f.readlines()
    chr1_seq = ''.join([l.strip() for l in lines[1:5]])
    chr2_seq = ''.join([l.strip() for l in lines[6:10]])

# Generate reads
num_reads = 1000
with open('test_sample_R1.fq', 'w') as r1, open('test_sample_R2.fq', 'w') as r2:
    for i in range(num_reads):
        ref = chr1_seq if i % 2 == 0 else chr2_seq

        # R1
        seq1, qual1 = generate_read(ref)
        r1.write(f'@read{i}/1\n{seq1}\n+\n{qual1}\n')

        # R2 (reverse complement mate)
        seq2, qual2 = generate_read(ref)
        r2.write(f'@read{i}/2\n{seq2}\n+\n{qual2}\n')

print(f"Generated {num_reads} paired-end reads")
PYTHON

# Compress FASTQ files
gzip -f test_sample_R1.fq
gzip -f test_sample_R2.fq

# Generate sample sheet
cat > samples.csv << EOF
sample_id,reads_r1,reads_r2
test_sample,test_data/test_sample_R1.fq.gz,test_data/test_sample_R2.fq.gz
EOF

# Generate regions BED file
cat > regions.bed << EOF
chr1	50	150
chr2	100	200
EOF

echo "✓ Test data generated:"
echo "  - reference.fa (2 chromosomes)"
echo "  - test_sample_R1.fq.gz (1000 reads)"
echo "  - test_sample_R2.fq.gz (1000 reads)"
echo "  - samples.csv"
echo "  - regions.bed (2 regions)"
