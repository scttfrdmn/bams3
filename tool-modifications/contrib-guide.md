# Contributing S3 Support Upstream

## Overview

This guide helps you contribute S3 support improvements to open source genomics tools.

## Target Projects

### Priority 1: htslib
**Why:** Foundation for entire ecosystem
**Repository:** https://github.com/samtools/htslib
**Status:** Has basic S3 support, needs improvements
**Impact:** HIGH - affects samtools, bcftools, and many tools

**Contribution opportunities:**
- [ ] Improve S3 plugin performance (prefetching, larger buffers)
- [ ] Better error messages and retry logic
- [ ] Write support improvements
- [ ] Documentation for end users

### Priority 2: samtools/bcftools
**Why:** Most widely used tools
**Repository:** https://github.com/samtools/samtools
**Status:** Uses htslib, mostly works
**Impact:** HIGH - direct user benefit

**Contribution opportunities:**
- [ ] Add S3-specific optimizations
- [ ] Document S3 usage patterns
- [ ] Add S3 URI examples to man pages

### Priority 3: BWA
**Why:** Standard read aligner
**Repository:** https://github.com/lh3/bwa
**Status:** No S3 support
**Impact:** MEDIUM - widely used but alternatives exist

**Contribution opportunities:**
- [ ] Add S3 support for reference genome
- [ ] Add S3 support for FASTQ input
- [ ] Maintain backward compatibility

### Priority 4: minimap2
**Why:** Modern aligner, gaining adoption
**Repository:** https://github.com/lh3/minimap2
**Status:** No S3 support
**Impact:** MEDIUM - growing user base

**Contribution opportunities:**
- [ ] Similar to BWA
- [ ] Author (Heng Li) is same as BWA

## Contribution Process

### 1. Understand the Project

**Before writing code:**

```bash
# Clone and explore
git clone https://github.com/samtools/htslib
cd htslib

# Read contribution guidelines
cat CONTRIBUTING.md

# Build and test
autoreconf -i
./configure --enable-s3
make
make test

# Understand existing S3 code
ls hfile_s3*
grep -r "s3://" .
```

**Questions to answer:**
- Who maintains the project?
- What's the coding style?
- Are there tests? How to run them?
- What's the release cycle?
- Are there open issues related to S3?

### 2. Discuss First

**Don't surprise maintainers with large PRs.**

**Good approach:**
1. Open an issue describing the problem
2. Propose your solution approach
3. Get feedback before implementing
4. Reference the issue in your PR

**Example issue:**

```markdown
Title: Improve S3 read performance with prefetching

Description:
I've been using htslib with S3 and noticed that sequential reads
are much slower than they could be. Each read() call makes a
separate S3 request.

I'd like to add prefetching that:
- Fetches larger chunks (e.g., 1MB) instead of the exact requested size
- Buffers ahead for sequential access patterns
- Falls back to current behavior for random access

This would improve performance for common use cases like:
- `samtools view` on entire file
- `bcftools query` processing

Would the maintainers be interested in such a contribution?
I'm happy to provide benchmarks and write tests.

Thoughts on the best approach?
```

### 3. Make Small, Focused Changes

**Bad PR:**
- 5000 lines changed
- Multiple features
- Refactoring + new features
- No tests

**Good PR:**
- <500 lines changed
- One clear improvement
- Tests included
- Documentation updated

**Break large work into multiple PRs:**

```
PR #1: Add prefetch buffer to hfile_s3.c
PR #2: Tune buffer sizes for S3
PR #3: Add retry logic for failed requests
PR #4: Update documentation
```

### 4. Follow Project Standards

**Code style:**
```c
// Check existing code style
// Most C projects use K&R or similar

// Example from htslib:
int hread(hFILE *fp, void *buffer, size_t nbytes)
{
    // 4-space indents
    // Opening brace on same line
    if (fp == NULL) {
        return -1;
    }
}
```

**Commit messages:**
```bash
# Good commit message format:
# <component>: <short description>
#
# Longer explanation if needed.
# - Bullet points OK
# - Reference issues with #123

# Example:
git commit -m "hfile_s3: Add prefetching for sequential reads

Improves performance for sequential S3 access by fetching larger
chunks (1MB) and buffering ahead. Falls back to current behavior
for random access.

Benchmarks show 3x improvement for 'samtools view' on S3 BAM files.

Fixes #12345"
```

### 5. Write Tests

**Every PR should include tests.**

```c
// Example test for S3 functionality
// test/s3_test.c

#include "htslib/hfile.h"
#include "tap.h"  // Test framework

int main(int argc, char **argv)
{
    plan_tests(5);

    // Test 1: Open S3 file
    hFILE *fp = hopen("s3://test-bucket/test.txt", "r");
    ok(fp != NULL, "Open S3 file");

    // Test 2: Read data
    char buf[100];
    ssize_t n = hread(fp, buf, sizeof(buf));
    ok(n > 0, "Read from S3 file");

    // Test 3: Seek (if supported)
    ok(hseek(fp, 10, SEEK_SET) == 10, "Seek in S3 file");

    // Test 4: Read after seek
    n = hread(fp, buf, 10);
    ok(n == 10, "Read after seek");

    // Test 5: Close
    ok(hclose(fp) == 0, "Close S3 file");

    return exit_status();
}
```

### 6. Benchmark and Document

**Include benchmarks in PR description:**

```markdown
## Performance Impact

Tested with 10GB BAM file on S3:

| Operation | Before | After | Improvement |
|-----------|--------|-------|-------------|
| Sequential read | 120s | 40s | 3x faster |
| Random access | 2s | 2s | No change |
| Memory usage | 64KB | 1.5MB | Acceptable |

Test setup:
- EC2 c5.4xlarge in us-east-1
- S3 bucket in us-east-1
- File: 1000genomes BAM (10GB)

Command: `time samtools view s3://bucket/file.bam > /dev/null`
```

### 7. Be Responsive

**After submitting PR:**
- Respond to comments within 24-48 hours
- Be open to changes
- Don't be defensive
- Explain your reasoning but accept feedback
- Make requested changes promptly

**Example response:**
```markdown
> Consider making the buffer size configurable

Good idea! I'll add an environment variable S3_BUFFER_SIZE
with a sensible default of 1MB. Will push an update today.
```

## Real Example: htslib S3 Prefetch

### Issue

Sequential reads of S3 files are slow due to many small requests.

### Solution Approach

1. Create GitHub issue (discuss with maintainers)
2. Implement prefetch buffer
3. Add tests
4. Benchmark
5. Submit PR

### Code Outline

```c
// hfile_s3_write.c additions

typedef struct {
    // Existing S3 handle fields...

    // New prefetch buffer
    unsigned char *prefetch_buf;
    size_t prefetch_size;      // Size of buffer (1MB default)
    size_t prefetch_pos;       // Current position in buffer
    size_t prefetch_fill;      // How much data in buffer
    off_t prefetch_offset;     // File offset of buffered data
} hFILE_s3;

static ssize_t s3_read(hFILE *fpv, void *buffer, size_t nbytes)
{
    hFILE_s3 *fp = (hFILE_s3 *)fpv;

    // Check if request can be satisfied from prefetch buffer
    if (fp->prefetch_buf &&
        fp->at >= fp->prefetch_offset &&
        fp->at < fp->prefetch_offset + fp->prefetch_fill) {

        // Serve from buffer
        size_t available = fp->prefetch_fill -
                          (fp->at - fp->prefetch_offset);
        size_t to_copy = nbytes < available ? nbytes : available;

        memcpy(buffer, fp->prefetch_buf + (fp->at - fp->prefetch_offset),
               to_copy);

        fp->at += to_copy;
        return to_copy;
    }

    // Buffer miss - fetch new chunk
    // Fetch larger chunk than requested for prefetching
    size_t fetch_size = nbytes > fp->prefetch_size ?
                        nbytes : fp->prefetch_size;

    // Make S3 request for larger chunk...
    // Store in prefetch buffer...
    // Return requested portion...
}
```

### PR Template

```markdown
# S3 prefetch buffer for improved sequential read performance

## Summary
Adds prefetching to S3 reads to reduce the number of requests for
sequential access patterns.

## Changes
- Added prefetch buffer (default 1MB) to hFILE_s3 structure
- Modified s3_read() to check buffer before making S3 requests
- Buffer size configurable via S3_BUFFER_SIZE env var
- Falls back to direct requests for random access

## Testing
- Added test_s3_prefetch to test suite
- Verified backward compatibility with existing tests
- All tests pass

## Performance
[Include benchmark table as shown above]

## Documentation
- Updated hfile.h with new behavior notes
- Updated S3 section in manual

## Checklist
- [x] Tests pass
- [x] Follows coding style
- [x] Documentation updated
- [x] Backward compatible
- [x] Benchmarked

Fixes #12345
```

## Upstream Acceptance Tips

### Do's
- ✅ Discuss before implementing
- ✅ Make changes optional (compile flag/env var)
- ✅ Keep backward compatibility
- ✅ Include tests
- ✅ Follow existing code style
- ✅ Write clear commit messages
- ✅ Be patient and responsive

### Don'ts
- ❌ Submit large PRs without discussion
- ❌ Break existing functionality
- ❌ Require new dependencies without discussion
- ❌ Ignore code review feedback
- ❌ Mix multiple features in one PR
- ❌ Skip tests

## Alternative: Maintain Fork

If upstream isn't interested or moves too slowly:

```bash
# Create your own maintained fork
git clone https://github.com/samtools/htslib
cd htslib
git remote add myorg https://github.com/myorg/htslib

# Make improvements
git checkout -b s3-improvements
# ... make changes ...
git commit -m "S3 improvements"
git push myorg s3-improvements

# Create releases
git tag v1.19.1-s3-optimized
git push myorg v1.19.1-s3-optimized

# Share with community
# Document improvements
# Provide binaries
```

**Document clearly:**
- What's different from upstream
- How to build and use
- How to contribute back
- When you'll merge upstream changes

## Success Stories

### Real examples of accepted contributions:
- htslib libcurl plugin (HTTP/HTTPS support)
- samtools CRAM support
- bcftools various format improvements

**Common pattern:**
1. User has need
2. Implements locally
3. Discusses with maintainers
4. Refines implementation
5. PR accepted
6. Everyone benefits

## Resources

- htslib contributing guide: https://github.com/samtools/htslib/blob/develop/CONTRIBUTING.md
- samtools mailing list: https://lists.sourceforge.net/lists/listinfo/samtools-help
- Code review checklist: https://github.com/samtools/htslib/wiki

## Questions?

Open an issue in this repository and we'll help you navigate the contribution process.
