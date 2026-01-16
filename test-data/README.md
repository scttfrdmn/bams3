# BAMS3 Test Data

Example data and test scripts for BAMS3 format validation.

## Quick Start

```bash
# From the aws-direct-s3 root directory:

# 1. Setup environment (using uv)
uv venv
uv pip install pysam boto3 pyarrow

# 2. Generate test data
cd test-data
../.venv/bin/python create_test_bam.py

# 3. Convert to BAMS3
../.venv/bin/python ../format-tools/bams3/bams3_converter.py \
    test_sample.bam test_sample.bams3

# 4. Query BAMS3
../.venv/bin/python ../format-tools/bams3/bams3_query.py \
    test_sample.bams3 chr1:1000000-2000000 --head 10

# 5. Get statistics
../.venv/bin/python ../format-tools/bams3/bams3_query.py \
    test_sample.bams3 --stats
```

## Test with S3

```bash
# Set your bucket (us-west-2 region)
export BAMS3_TEST_BUCKET=your-bucket-name

# Run full S3 integration test
./test_s3_workflow.sh
```

## Files

- `create_test_bam.py` - Generate synthetic BAM files for testing
- `test_s3_workflow.sh` - Complete S3 integration test
- `TEST_RESULTS.md` - Detailed test results and performance analysis
- `test_sample.bam` - Generated test BAM (created by script)
- `test_sample.bams3/` - Generated BAMS3 dataset (created by script)

## Test Results

See [TEST_RESULTS.md](TEST_RESULTS.md) for complete test results.

**Summary:**
- ✅ All tools working
- ✅ 7-50,000x less data transfer for region queries
- ✅ 18,000x faster statistics (metadata vs scan)
- ✅ Ready for production Go/Rust implementation

## Creating Larger Test Data

```bash
# Create larger test file (10,000 reads)
../.venv/bin/python create_test_bam.py large_test.bam 10000

# Convert with custom chunk size (5Mbp)
../.venv/bin/python ../format-tools/bams3/bams3_converter.py \
    large_test.bam large_test.bams3 --chunk-size 5000000
```

## Benchmarking

Compare BAMS3 vs traditional BAM access:

```bash
# Query with traditional BAM
time samtools view test_sample.bam chr1:1000000-2000000 > /dev/null

# Query with BAMS3
time ../.venv/bin/python ../format-tools/bams3/bams3_query.py \
    test_sample.bams3 chr1:1000000-2000000 > /dev/null

# For small files, difference is minimal
# For large files (GB+), BAMS3 is orders of magnitude faster
```

## AWS Configuration

The test scripts use:
- **Profile**: `aws` (default AWS CLI profile)
- **Region**: `us-west-2`

If you need different settings:

```bash
# Use different profile
export AWS_PROFILE=your-profile

# Or set in script
AWS_PROFILE=your-profile ./test_s3_workflow.sh
```

## Next Steps

1. ✅ Test Python POC (done - see TEST_RESULTS.md)
2. 🔧 Implement Go version (see ../format-tools/bams3/GO_RUST_IMPLEMENTATIONS.md)
3. 🔧 Add binary chunk format (replace JSON)
4. 🔧 Add compression (zstd/lz4)
5. 🚀 Production deployment
