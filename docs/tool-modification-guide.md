# Tool Modification Guide: Adding S3 Support

## Overview

This guide covers practical approaches for adding S3 support to existing bioinformatics tools. The right approach depends on the tool's language, architecture, and access patterns.

## Decision Tree

```
Is the tool written in C/C++?
├─ Yes: Does it use htslib?
│  ├─ Yes → Use htslib's S3 support (easiest)
│  └─ No → Add AWS SDK or libcurl
└─ No: What language?
   ├─ Python → Use boto3/s3fs
   ├─ Java → Use AWS SDK for Java
   ├─ R → Use aws.s3 package
   └─ Other → Check for AWS SDK support
```

## Approach 1: Leverage Existing Libraries

### For htslib-based tools (C/C++)

**Best case:** Tool already uses htslib for I/O.

**Changes needed:** NONE if htslib built with S3 support!

**Example:**
```c
// This code already works with S3 if htslib has S3 enabled:
htsFile *fp = hts_open("s3://bucket/file.bam", "r");
// All existing htslib calls work unchanged
```

**Action items:**
1. Verify tool uses htslib for file I/O
2. Ensure htslib built with `--enable-s3`
3. Update documentation to mention S3 URIs work
4. Test with S3 URLs

**Tools in this category:**
- samtools, bcftools (already work!)
- Any tool using htslib API

### For Python tools

**Use existing S3-aware libraries:**

```python
# Option 1: boto3 (low-level)
import boto3
s3 = boto3.client('s3')
obj = s3.get_object(Bucket='bucket', Key='file.txt')
data = obj['Body'].read()

# Option 2: s3fs (filesystem interface)
import s3fs
fs = s3fs.S3FileSystem()
with fs.open('s3://bucket/file.txt') as f:
    data = f.read()

# Option 3: smart_open (handles s3:// URLs transparently)
from smart_open import open
with open('s3://bucket/file.txt', 'rb') as f:
    data = f.read()
```

**pandas/dask already support S3:**
```python
import pandas as pd
df = pd.read_csv('s3://bucket/file.csv')
df = pd.read_parquet('s3://bucket/file.parquet')
```

### For Java tools

**Use AWS SDK for Java:**

```java
import software.amazon.awssdk.services.s3.S3Client;
import software.amazon.awssdk.services.s3.model.GetObjectRequest;

S3Client s3 = S3Client.builder().build();
GetObjectRequest request = GetObjectRequest.builder()
    .bucket("bucket")
    .key("file.txt")
    .build();
InputStream stream = s3.getObject(request);
```

**Note:** GATK already has some Google Cloud support; similar patterns could work for AWS.

## Approach 2: Virtual Filesystem Layer

### When to use
- Tool has many file operations
- Refactoring all file I/O is impractical
- Want minimal code changes

### Implementation pattern

Create an abstraction layer:

```c
// vfs.h - Virtual filesystem interface
typedef struct {
    void* handle;
    int (*read)(void* handle, char* buf, size_t size);
    int (*write)(void* handle, const char* buf, size_t size);
    int (*seek)(void* handle, off_t offset, int whence);
    int (*close)(void* handle);
} vfs_file;

vfs_file* vfs_open(const char* uri, const char* mode);

// vfs.c - Implementation
vfs_file* vfs_open(const char* uri, const char* mode) {
    if (strncmp(uri, "s3://", 5) == 0) {
        return s3_open(uri, mode);
    } else {
        return posix_open(uri, mode);
    }
}
```

**Then replace throughout codebase:**
```c
// Before:
FILE* fp = fopen(path, "r");

// After:
vfs_file* fp = vfs_open(path, "r");
```

### Example: BWA modification

BWA reads FASTA reference and FASTQ inputs. Key functions:
- `err_fopen()` - Opens files
- `fread()` - Reads data
- `bwa_idx_load()` - Loads reference index

**Modification approach:**
1. Create VFS layer (above)
2. Replace file operations in key functions
3. Add S3 backend using AWS SDK or libcurl
4. Maintain backward compatibility

## Approach 3: Streaming/Piping

### When to use
- Tool reads input sequentially
- Tool accepts stdin input
- Don't want to modify tool source

### Implementation

Use external S3 streaming:

```bash
# Stream FASTQ from S3 to tool stdin
aws s3 cp s3://bucket/file.fastq.gz - | gunzip | tool -

# Or with custom streaming tool:
s3-stream s3://bucket/file.fastq.gz | tool -
```

**Create a simple s3-stream utility:**

```python
#!/usr/bin/env python3
import sys
import boto3

def stream_s3_to_stdout(s3_uri):
    # Parse s3://bucket/key
    parts = s3_uri[5:].split('/', 1)
    bucket, key = parts[0], parts[1]

    s3 = boto3.client('s3')
    obj = s3.get_object(Bucket=bucket, Key=key)

    # Stream to stdout
    for chunk in obj['Body'].iter_chunks():
        sys.stdout.buffer.write(chunk)

if __name__ == '__main__':
    stream_s3_to_stdout(sys.argv[1])
```

**Limitations:**
- Only works for sequential read
- No random access
- No write support

## Approach 4: Build-time Configuration

### For tools with plugin architecture

Some tools support plugins or can be extended:

```c
// Tool provides plugin hooks
struct io_plugin {
    const char* scheme;  // "s3"
    void* (*open)(const char* uri, const char* mode);
    int (*read)(void* handle, void* buf, size_t size);
    int (*close)(void* handle);
};

// Register S3 plugin at runtime
register_io_plugin(&s3_plugin);
```

**Advantages:**
- No changes to core tool
- Plugins can be optional
- Easy to maintain separately

**Example:** htslib's approach!

## Implementation Checklist

### Before Starting
- [ ] Identify all file I/O operations in tool
- [ ] Check if tool already uses abstraction library (htslib, etc.)
- [ ] Determine access pattern (sequential, random, read, write)
- [ ] Check for existing plugin/extension mechanism

### Core Functionality
- [ ] URI parsing (s3://bucket/key)
- [ ] Credential handling (AWS SDK default chain)
- [ ] Read support with streaming
- [ ] Range request support (for random access)
- [ ] Error handling and retries

### Optimization
- [ ] Buffer size tuning (MB not KB for S3)
- [ ] Prefetching for sequential reads
- [ ] Parallel fetching for large files
- [ ] Caching for repeated access

### Testing
- [ ] Unit tests for S3 operations
- [ ] Integration tests with real S3
- [ ] Performance benchmarks vs local files
- [ ] Error condition handling

### Documentation
- [ ] How to configure credentials
- [ ] URI format examples
- [ ] Performance characteristics
- [ ] Known limitations

### Upstream Contribution
- [ ] Follow project coding style
- [ ] Make S3 support optional (compile flag)
- [ ] Don't break existing functionality
- [ ] Include tests
- [ ] Update documentation

## Common Pitfalls

### 1. Small Buffer Sizes
```c
// Bad: Too small for S3
char buf[4096];

// Good: Larger buffers for cloud
char buf[1024 * 1024];  // 1MB
```

### 2. No Retry Logic
```c
// Bad: One attempt
response = s3_get_object(bucket, key);

// Good: Retry with exponential backoff
for (int retry = 0; retry < 3; retry++) {
    response = s3_get_object(bucket, key);
    if (response.success) break;
    sleep(pow(2, retry));  // 1s, 2s, 4s
}
```

### 3. Synchronous Large Downloads
```c
// Bad: Download entire file first
download_file(s3_uri, "/tmp/file");
process_file("/tmp/file");

// Good: Stream as you process
stream = s3_open(s3_uri);
while (chunk = stream_read(stream, buf, sizeof(buf))) {
    process_chunk(buf, chunk);
}
```

### 4. Not Using Range Requests
```c
// Bad: Download entire BAM to read one region
download_entire_bam(s3_uri);

// Good: Use range requests with index
bam_index = read_bai(s3_uri + ".bai");
chunk_offset = index_lookup(bam_index, region);
chunk = s3_get_range(s3_uri, chunk_offset, chunk_size);
```

## Language-Specific Resources

### C/C++
- AWS SDK for C++: https://github.com/aws/aws-sdk-cpp
- libcurl for simple HTTP: https://curl.se/libcurl/
- htslib S3 plugin: Reference implementation

### Python
- boto3: https://boto3.amazonaws.com/
- s3fs: https://s3fs.readthedocs.io/
- smart_open: https://github.com/RaRe-Technologies/smart_open

### Java
- AWS SDK for Java: https://aws.amazon.com/sdk-for-java/

### R
- aws.s3: https://github.com/cloudyr/aws.s3
- Arrow with S3: https://arrow.apache.org/docs/r/

## Next Steps

1. Choose target tool for modification
2. Analyze its I/O patterns
3. Select appropriate approach from above
4. Create prototype
5. Benchmark performance
6. Refine implementation
7. Contribute upstream

See `tool-modifications/` directory for specific examples.
