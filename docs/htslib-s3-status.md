# HTSlib S3 Support - Current Status

## Overview

HTSlib is the foundational C library for reading/writing genomics file formats (SAM/BAM/CRAM/VCF/BCF). It's used by:
- samtools
- bcftools
- Many other tools in the ecosystem

**Key insight:** Adding S3 support to htslib propagates to the entire ecosystem.

## Current State

### S3 Support Exists via Plugins

HTSlib has a **plugin architecture** for remote file access:
- `hfile_s3.c` - S3 plugin using libcurl + AWS SDK
- `hfile_libcurl.c` - Generic HTTP/HTTPS support

**Compilation requirement:**
```bash
# Build with S3 support
./configure --enable-s3 --enable-libcurl
make
```

### What Works Today

```bash
# These already work if htslib built with S3 support:
samtools view s3://1000genomes/data/file.bam

# Set credentials via environment or AWS config
export AWS_ACCESS_KEY_ID=...
export AWS_SECRET_ACCESS_KEY=...
# Or use AWS SDK default credential chain
```

**Supported operations:**
- Read BAM/CRAM/VCF from S3
- Range requests for indexed formats (BAM + BAI)
- Random access when index available

## Known Limitations and Gaps

### 1. Performance Issues
- **No prefetching** - Each read is a separate S3 request
- **No parallel fetch** - Single-threaded S3 access
- **Small buffer sizes** - Optimized for local files, not cloud
- **No intelligent caching** - Doesn't learn from access patterns

### 2. Write Support Limited
- Reading is well-supported
- **Writing to S3 is incomplete/buggy**
- Multipart upload not properly implemented
- No atomic writes or error recovery

### 3. Configuration Complexity
- Not built by default in many distributions
- Users must compile from source with flags
- Credentials setup non-obvious
- No clear documentation for end users

### 4. Error Handling
- Network errors not always handled gracefully
- Retry logic basic or missing
- Timeout handling needs improvement

### 5. Missing Features
- No support for S3 Select (server-side filtering)
- No integration with AWS S3 Transfer Acceleration
- No support for S3 event notifications
- Limited support for S3-compatible endpoints (MinIO, Ceph)

## Opportunities for Improvement

### High Impact
1. **Optimize for S3 access patterns**
   - Larger read buffers (MB not KB)
   - Prefetch subsequent blocks
   - Parallel range requests
   - Adaptive caching based on access patterns

2. **Improve write support**
   - Proper multipart upload implementation
   - Atomic writes with abort on error
   - Progress reporting for large files

3. **Better error handling**
   - Exponential backoff retries
   - Handle throttling (429) gracefully
   - Detailed error messages for users

### Medium Impact
4. **Distribution improvements**
   - Compile S3 support by default
   - Better documentation
   - Pre-built binaries with S3 enabled

5. **Advanced features**
   - S3 Select integration for VCF queries
   - Support for S3 Glacier retrieval
   - CloudFront integration for public datasets

### Research Needed
- Benchmark current S3 performance vs local file
- Profile to find bottlenecks
- Compare with other tools (e.g., Google Cloud's implementation)

## Related Work

### Other Cloud Storage Plugins
- **Google Cloud Storage** - `hfile_gcs.c`
  - More mature implementation
  - Good reference for S3 improvements

- **Azure Blob Storage** - Limited/no native support

### Alternative Approaches
- **FUSE layer** - Mount S3, htslib sees local files
  - Works but adds overhead
  - Doesn't optimize for genomics access patterns

- **Wrapper libraries** - Intercept file calls
  - Complex to maintain
  - Can't optimize at application level

## Action Items for This Project

### Phase 1: Assessment
- [ ] Build htslib with S3 support enabled
- [ ] Test current S3 functionality with real data
- [ ] Benchmark read performance (local vs S3)
- [ ] Document current behavior and limitations

### Phase 2: Improvements
- [ ] Implement prefetching for sequential access
- [ ] Optimize buffer sizes for cloud
- [ ] Add better error messages
- [ ] Improve retry logic

### Phase 3: New Features
- [ ] Parallel range fetching
- [ ] Write support improvements
- [ ] S3-specific optimizations

### Phase 4: Upstream Contribution
- [ ] Clean up patches
- [ ] Write tests
- [ ] Submit PRs to htslib project
- [ ] Work with maintainers on integration

## Resources

- HTSlib GitHub: https://github.com/samtools/htslib
- HTSlib plugin architecture docs
- AWS SDK for C++
- HTSJDK (Java implementation) - has good cloud support for reference

## Notes

The htslib maintainers are receptive to improvements but careful about dependencies and portability. Any contributions should:
- Be optional (compile-time flag)
- Not break existing functionality
- Work across platforms
- Include tests
- Follow coding style
