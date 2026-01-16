# FUSE Mounting Examples

## Overview

FUSE (Filesystem in Userspace) lets you mount S3 as a local filesystem. Tools see regular POSIX paths and work unchanged.

**Pros:**
- Zero tool modifications needed
- Works with any POSIX application
- Easy for naive users

**Cons:**
- Performance overhead
- Caching behavior critical
- Can be unpredictable

## Options Compared

| Tool | Speed | Maturity | Features | Best For |
|------|-------|----------|----------|----------|
| s3fs-fuse | ⭐⭐ | High | Full POSIX | Compatibility |
| goofys | ⭐⭐⭐⭐ | Medium | Limited writes | Read-heavy |
| mountpoint-for-amazon-s3 | ⭐⭐⭐⭐⭐ | New (AWS) | Read-only | Best performance |

## Option 1: mountpoint-for-amazon-s3 (Recommended)

AWS's official solution, written in Rust, optimized for high-throughput read access.

### Installation

```bash
# macOS
brew install mountpoint-s3

# Linux (Amazon Linux 2023)
sudo yum install mountpoint-s3

# Linux (Ubuntu/Debian)
wget https://s3.amazonaws.com/mountpoint-s3-release/latest/x86_64/mount-s3.deb
sudo apt install ./mount-s3.deb

# Or build from source
git clone https://github.com/awslabs/mountpoint-s3.git
cd mountpoint-s3
cargo build --release
```

### Basic Usage

```bash
# Create mount point
mkdir ~/s3-data

# Mount a bucket
mount-s3 my-genomics-bucket ~/s3-data

# Use it like local filesystem
ls ~/s3-data
samtools view ~/s3-data/sample1.bam | head

# Unmount
umount ~/s3-data
```

### Advanced Configuration

```bash
# High-performance options
mount-s3 \
    --cache /tmp/s3-cache \            # Local cache directory
    --max-cache-size 100000 \          # 100GB cache
    --allow-other \                    # All users can access
    --read-only \                      # Prevent accidental writes
    --prefix data/genomics/ \          # Mount subdirectory only
    my-bucket \
    ~/s3-data

# For parallel processing (multiple processes)
mount-s3 \
    --cache /tmp/s3-cache \
    --thread-count 16 \                # More threads for parallel I/O
    --part-size 8388608 \              # 8MB parts
    my-bucket \
    ~/s3-data
```

### Credentials

```bash
# Option 1: AWS CLI configured
aws configure

# Option 2: Environment variables
export AWS_ACCESS_KEY_ID=...
export AWS_SECRET_ACCESS_KEY=...
export AWS_DEFAULT_REGION=us-east-1

# Option 3: IAM role (on EC2)
# No configuration needed - uses instance role
```

### Limitations

- **Read-only** - Cannot write files (by design)
- **No listing caching** - `ls` is slower than local
- **Eventual consistency** - May not see recent writes immediately

### Performance Tips

1. **Enable caching:**
   ```bash
   mount-s3 --cache /fast-nvme/cache --max-cache-size 500000 bucket mount/
   ```

2. **Increase part size for large files:**
   ```bash
   mount-s3 --part-size 16777216 bucket mount/  # 16MB parts
   ```

3. **Use more threads:**
   ```bash
   mount-s3 --thread-count $(nproc) bucket mount/
   ```

## Option 2: goofys

Fast S3 filesystem optimized for read-heavy workloads.

### Installation

```bash
# macOS
brew install goofys

# Linux
wget https://github.com/kahing/goofys/releases/latest/download/goofys
chmod +x goofys
sudo mv goofys /usr/local/bin/
```

### Usage

```bash
# Basic mount
goofys my-bucket ~/s3-data

# With caching
goofys \
    --stat-cache-ttl 60s \      # Cache stat() calls
    --type-cache-ttl 60s \      # Cache type lookups
    --dir-mode 0755 \
    --file-mode 0644 \
    my-bucket \
    ~/s3-data
```

### Best For

- Read-heavy workloads
- Need better performance than s3fs
- Can tolerate limited write support

## Option 3: s3fs-fuse

Most mature option, full POSIX support, slower but stable.

### Installation

```bash
# macOS
brew install s3fs

# Ubuntu/Debian
sudo apt-get install s3fs

# Build from source
git clone https://github.com/s3fs-fuse/s3fs-fuse.git
cd s3fs-fuse
./autogen.sh
./configure
make
sudo make install
```

### Usage

```bash
# Create credentials file
echo ACCESS_KEY_ID:SECRET_ACCESS_KEY > ~/.passwd-s3fs
chmod 600 ~/.passwd-s3fs

# Mount
s3fs my-bucket ~/s3-data -o passwd_file=~/.passwd-s3fs

# With caching and options
s3fs my-bucket ~/s3-data \
    -o passwd_file=~/.passwd-s3fs \
    -o use_cache=/tmp/s3fs-cache \
    -o max_stat_cache_size=10000 \
    -o stat_cache_expire=60 \
    -o multireq_max=500 \
    -o parallel_count=20
```

### Best For

- Need write support
- Need full POSIX semantics
- Performance is not critical

## Genomics Pipeline Example

### Pipeline with mountpoint-s3

```bash
#!/bin/bash
# genomics-pipeline-fuse.sh

# Setup
BUCKET="my-genomics-bucket"
MOUNT_POINT="$HOME/s3-mount"
CACHE_DIR="/mnt/nvme/s3-cache"

# Mount S3
mkdir -p $MOUNT_POINT $CACHE_DIR
mount-s3 $BUCKET $MOUNT_POINT \
    --cache $CACHE_DIR \
    --max-cache-size 200000 \
    --thread-count 16

# Reference genome (frequently accessed - will be cached)
REF="$MOUNT_POINT/references/hg38.fa"

# Input FASTQ files
R1="$MOUNT_POINT/fastq/sample1_R1.fastq.gz"
R2="$MOUNT_POINT/fastq/sample1_R2.fastq.gz"

# Output (write to local disk, then upload)
OUTPUT_DIR="/tmp/analysis"
mkdir -p $OUTPUT_DIR

# Run pipeline
echo "Running BWA alignment..."
bwa mem -t 16 $REF $R1 $R2 | \
    samtools view -b -o $OUTPUT_DIR/aligned.bam

echo "Sorting BAM..."
samtools sort -@ 16 $OUTPUT_DIR/aligned.bam -o $OUTPUT_DIR/sorted.bam

echo "Indexing BAM..."
samtools index $OUTPUT_DIR/sorted.bam

# Upload results
echo "Uploading results..."
aws s3 cp $OUTPUT_DIR/sorted.bam s3://$BUCKET/results/
aws s3 cp $OUTPUT_DIR/sorted.bam.bai s3://$BUCKET/results/

# Cleanup
rm -rf $OUTPUT_DIR
umount $MOUNT_POINT
```

### Auto-mount on Boot (systemd)

```ini
# /etc/systemd/system/s3-mount.service
[Unit]
Description=Mount S3 bucket for genomics
After=network.target

[Service]
Type=forking
User=genomics
ExecStart=/usr/local/bin/mount-s3 \
    --cache /mnt/cache/s3 \
    --max-cache-size 500000 \
    --allow-other \
    my-genomics-bucket \
    /data/s3
ExecStop=/bin/umount /data/s3
Restart=on-failure

[Install]
WantedBy=multi-user.target
```

```bash
# Enable and start
sudo systemctl enable s3-mount
sudo systemctl start s3-mount

# Check status
sudo systemctl status s3-mount
```

## Benchmarking

### Test Script

```bash
#!/bin/bash
# benchmark-fuse.sh

MOUNT_POINT="$HOME/s3-mount"
TEST_FILE="$MOUNT_POINT/test-data/large-file.bam"  # 10GB file

echo "Testing FUSE mount performance..."

# Test 1: Sequential read
echo "Sequential read:"
time samtools view $TEST_FILE > /dev/null

# Test 2: Random access by region
echo "Random access (with index):"
time samtools view $TEST_FILE chr1:1000000-2000000 > /dev/null

# Test 3: Multiple small accesses
echo "Multiple small accesses:"
time for i in {1..10}; do
    samtools view $TEST_FILE chr1:${i}000000-${i}100000 > /dev/null
done

# Compare with direct S3 access (if tool supports it)
echo "Direct S3 access:"
time samtools view s3://my-bucket/test-data/large-file.bam > /dev/null
```

### Expected Results

| Operation | Local SSD | mountpoint-s3 | goofys | s3fs |
|-----------|-----------|---------------|--------|------|
| Sequential read (10GB) | 10s | 30s | 60s | 120s |
| Random access (1MB) | 0.1s | 0.5s | 1s | 2s |
| Stat file | <0.01s | 0.1s | 0.1s | 0.2s |
| List directory (1000 files) | 0.05s | 2s | 1s | 5s |

## Troubleshooting

### Mount fails

```bash
# Check credentials
aws s3 ls s3://my-bucket/

# Check mount point exists and is empty
mkdir -p ~/s3-mount
ls ~/s3-mount  # Should be empty

# Try with debug output
mount-s3 --debug my-bucket ~/s3-mount
```

### Slow performance

```bash
# Check cache is being used
mount-s3 --cache /fast-disk/cache my-bucket ~/s3-mount

# Monitor cache usage
watch -n 1 du -sh /fast-disk/cache

# Increase threads
mount-s3 --thread-count 32 my-bucket ~/s3-mount
```

### "Transport endpoint not connected"

```bash
# Unmount forcefully
fusermount -u ~/s3-mount
# Or on macOS
umount -f ~/s3-mount

# Remount
mount-s3 my-bucket ~/s3-mount
```

## Next Steps

1. Install mountpoint-for-amazon-s3
2. Test with your data
3. Benchmark vs copying data locally
4. Tune cache settings for your workload
5. Consider if tool modification would be better long-term

See `../genomics-pipeline/` for complete workflow examples.
