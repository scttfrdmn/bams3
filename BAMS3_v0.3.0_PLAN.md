# BAMS3 v0.3.0 Implementation Plan

**Status:** Planning
**Target Date:** 2026-02-01 (2 weeks)
**Focus:** Performance, S3 Integration, Production Features

## Executive Summary

v0.2.0 achieved production-ready binary format with excellent query performance (92x faster) and storage efficiency (10% smaller than BAM). However, conversion speed (37,700 reads/sec) is slower than JSON format (51,600 reads/sec) due to complex binary encoding.

**v0.3.0 Goals:**
1. **3-5x faster conversion** through parallel compression
2. **Direct S3 access** for true cloud-native queries
3. **Streaming processing** for memory efficiency
4. **Production hardening** for large-scale deployments

## Performance Baseline (v0.2.0)

| Metric | Current | Target v0.3 | Improvement |
|--------|---------|-------------|-------------|
| Conversion | 37,700 reads/s | 100,000+ reads/s | 2.7x |
| Query (100KB) | 1.3s (local) | 2.0s (S3 direct) | 1.5x slower but no download |
| Memory usage | ~500 MB | ~100 MB | 5x reduction |
| S3 access | Download first | Stream directly | Eliminates download |

## Priority Features

### P0: Critical for Production

#### 1. Parallel Compression ⭐⭐⭐
**Goal:** 3-5x faster conversion through multi-core processing

**Current bottleneck:**
- Single-threaded conversion: 37,700 reads/sec
- Compression is CPU-intensive and parallelizable
- Not using all available cores (8+ cores idle)

**Implementation:**
1. **Parallel chunk compression** (`pkg/bams3/parallel_writer.go`)
   - Worker pool pattern with N goroutines (N = CPU cores)
   - Each worker compresses one chunk
   - Main thread coordinates chunk assignment

2. **Channel-based pipeline:**
   ```go
   reads (input) → chunking → compression workers → disk writers
   ```

3. **Memory management:**
   - Limit concurrent chunks to avoid OOM
   - Backpressure when memory limit reached
   - Reuse compression buffers

**Expected results:**
- 8-core machine: 8x speedup (theoretical)
- Practical: 3-5x speedup (due to I/O overhead)
- Target: 100,000+ reads/sec

**Effort:** 3-4 days
**Risk:** Low (well-understood pattern)

#### 2. Direct S3 Reading ⭐⭐⭐
**Goal:** Query BAMS3 datasets directly from S3 without local download

**Current limitation:**
- Must download entire dataset to local disk first
- Query: `aws s3 sync` (minutes) + `bams3 query` (seconds)
- Storage: Need local space for full dataset

**Implementation:**
1. **S3 reader** (`pkg/bams3/s3_reader.go`)
   - Use AWS SDK Go v2
   - Support S3 URI: `s3://bucket/dataset.bams3/`
   - Lazy loading: Fetch only queried chunks

2. **Metadata caching:**
   ```go
   // Cache metadata locally for fast subsequent queries
   ~/.cache/bams3/
   ├── datasets.db         # Dataset metadata cache
   └── chunks/             # Optional chunk cache
   ```

3. **Range requests:**
   - Use S3 byte-range requests if needed
   - Download only required chunk files
   - Stream decompression on-the-fly

**Expected results:**
- Query without download: `bams3 query s3://bucket/sample.bams3 chr1:1M-2M`
- First query: ~2-3s (vs current: minutes)
- Subsequent queries: ~1-2s (cached metadata)
- No local storage needed

**Effort:** 4-5 days
**Risk:** Medium (S3 error handling, retries, credentials)

#### 3. Streaming Query Processing ⭐⭐
**Goal:** Process reads without loading entire chunk into memory

**Current limitation:**
- Load entire chunk (10-40 MB) into memory
- Decompress all data before filtering
- Memory usage: 500+ MB for queries

**Implementation:**
1. **Streaming decompressor** (`pkg/bams3/stream_reader.go`)
   - Stream zstd decompression
   - Process records one-by-one
   - Discard non-matching reads immediately

2. **Memory-mapped reading:**
   - Use mmap for large chunks
   - Lazy decompression
   - Process in record-sized windows

**Expected results:**
- Memory usage: <100 MB (vs 500 MB)
- Slightly slower queries (~10% overhead)
- Support very large chunks (>100 MB)

**Effort:** 3-4 days
**Risk:** Low

### P1: Important for Production

#### 4. Chunk Metadata in Index ⭐
**Goal:** Skip chunks without reading them

**Implementation:**
- Add per-chunk statistics to spatial index:
  ```json
  {
    "path": "chr1/000000000-001048576.chunk",
    "stats": {
      "min_mapq": 0,
      "max_mapq": 60,
      "has_duplicates": false,
      "read_count": 50000
    }
  }
  ```

**Use cases:**
- Skip chunks with no high-quality reads (min_mapq filter)
- Skip chunks with no duplicates (dedup filter)
- Count reads without loading chunks

**Effort:** 2 days
**Risk:** Low

#### 5. Progress Reporting ⭐
**Goal:** Show progress during long-running operations

**Implementation:**
- Progress bar for conversion
- ETA for large files
- Cancellation support (Ctrl+C)
- Resume support for interrupted conversions

**Effort:** 2 days
**Risk:** Low

#### 6. Validation Tools ⭐
**Goal:** Verify dataset integrity

**Implementation:**
- `bams3 validate` command
- Check chunk checksums
- Verify metadata consistency
- Compare with original BAM

**Effort:** 2 days
**Risk:** Low

### P2: Nice to Have

#### 7. Configurable Compression Levels
**Options:** zstd levels 1-19
- Level 1: Fast, lower compression
- Level 3: Default (balanced)
- Level 19: Slow, maximum compression

**Effort:** 1 day

#### 8. Columnar Storage Option
**Experimental feature for coverage queries:**
- Store sequences separate from qualities
- Query only needed columns
- 10-100x less data for coverage queries

**Effort:** 1 week (experimental)

#### 9. Multi-Sample Datasets
**Support for population studies:**
- Multiple samples in one dataset
- Per-sample indexes
- Efficient multi-sample queries

**Effort:** 2 weeks (defer to v0.4)

## Implementation Timeline

### Week 1 (Jan 15-22)
- ✅ Day 1-2: Binary format validation (DONE)
- ✅ Day 3: v0.3.0 planning (DONE)
- 🔄 Day 4-5: Parallel compression implementation
- 📝 Day 6-7: Testing and benchmarking

### Week 2 (Jan 22-29)
- 📝 Day 8-10: Direct S3 reading implementation
- 📝 Day 11-12: Streaming query processing
- 📝 Day 13-14: Integration testing

### Week 3 (Jan 29 - Feb 1)
- 📝 Day 15-16: Performance optimization
- 📝 Day 17: Documentation and examples
- 📝 Day 18: Release v0.3.0

## Testing Strategy

### Unit Tests
- Parallel compression correctness
- S3 mocking for reader tests
- Streaming decompression validation
- Memory leak detection

### Integration Tests
- Real S3 bucket tests (with cleanup)
- Large file tests (>1 GB)
- Concurrent query tests
- Error recovery tests

### Performance Benchmarks
- Conversion throughput: 100K+ reads/sec target
- Query latency: <2s for S3 direct access
- Memory usage: <100 MB for queries
- Parallel scalability: Linear up to 8 cores

## Risk Mitigation

### High Risk: S3 Integration
**Risk:** Network errors, credential issues, S3 API limits
**Mitigation:**
- Comprehensive error handling
- Exponential backoff with retries
- Clear error messages
- Graceful degradation (fallback to local)

### Medium Risk: Parallel Compression
**Risk:** Race conditions, memory exhaustion, deadlocks
**Mitigation:**
- Thorough testing with race detector
- Memory limits and backpressure
- Cancellation and cleanup
- Worker pool best practices

### Low Risk: Streaming Processing
**Risk:** Correctness bugs, performance regression
**Mitigation:**
- Compare output with non-streaming
- Benchmark both paths
- Optional streaming flag initially

## Success Criteria

### Must Have (P0)
- ✅ Conversion: ≥100,000 reads/sec
- ✅ S3 queries work without local download
- ✅ Memory usage: ≤100 MB for queries
- ✅ All v0.2.0 features still work

### Should Have (P1)
- ✅ Progress reporting for conversions >1 minute
- ✅ Validation tools for integrity checks
- ✅ Chunk metadata in index

### Nice to Have (P2)
- ⚠️ Configurable compression levels
- ⚠️ Columnar storage (experimental)
- ❌ Multi-sample datasets (defer to v0.4)

## File Structure Changes

### New Files
```
bams3-go/
├── pkg/bams3/
│   ├── parallel_writer.go      # P0: Parallel compression
│   ├── s3_reader.go             # P0: S3 reading
│   ├── s3_writer.go             # P0: S3 writing
│   ├── stream_reader.go         # P0: Streaming queries
│   └── cache.go                 # P1: Metadata caching
├── cmd/bams3/
│   └── validate.go              # P1: Validation command
└── docs/
    ├── S3_INTEGRATION.md        # S3 usage guide
    └── PERFORMANCE_TUNING.md    # Optimization guide
```

### Modified Files
```
pkg/bams3/
├── writer.go          # Add parallel option
├── reader.go          # Add S3 support
└── types.go           # Add chunk stats
```

## Breaking Changes

### None Planned
v0.3.0 will be backward compatible with v0.2.0:
- All v0.2.0 datasets remain valid
- CLI flags unchanged (new flags added)
- Library API stable (new functions added)

## Documentation Updates

### User Documentation
1. S3 usage guide
   - AWS credentials configuration
   - S3 URI format
   - Caching behavior
   - Performance tips

2. Performance tuning
   - Chunk size selection
   - Compression level selection
   - Parallel workers configuration
   - Memory limits

3. Migration guide
   - Uploading to S3
   - Querying from S3
   - Cost optimization

### Developer Documentation
1. Architecture overview
2. Parallel processing patterns
3. S3 integration details
4. Testing guidelines

## Dependencies

### New Dependencies
- AWS SDK Go v2 (`github.com/aws/aws-sdk-go-v2`)
- AWS S3 module (`github.com/aws/aws-sdk-go-v2/service/s3`)
- AWS credentials module

### Updated Dependencies
- None (keep existing versions)

## Backward Compatibility

### v0.2.0 → v0.3.0
- ✅ All v0.2.0 datasets readable
- ✅ All v0.2.0 CLI commands work
- ✅ All v0.2.0 library APIs stable
- ➕ New features are additive only

### v0.1.0 → v0.3.0
- ✅ JSON format still supported (deprecated)
- ✅ Auto-detection still works
- ⚠️ Recommendation: Convert to binary format

## Release Checklist

- [ ] All P0 features implemented and tested
- [ ] Performance targets met (100K reads/sec)
- [ ] S3 integration working with real buckets
- [ ] Memory usage under target (100 MB)
- [ ] Documentation updated
- [ ] CHANGELOG.md updated
- [ ] Real data testing (chr22 + full WGS)
- [ ] EC2 testing (optional, validate latency)
- [ ] Version bumped to 0.3.0
- [ ] Git tag created
- [ ] Release notes published

## Post-Release

### Monitoring
- Community feedback on S3 performance
- Memory usage in production
- Conversion speed improvements
- Bug reports and fixes

### v0.4.0 Planning
Based on v0.3.0 learnings:
1. Columnar storage (if experimental successful)
2. Multi-sample datasets
3. Advanced indexing (Bloom filters)
4. Expression quantification support

## Questions for User

1. **EC2 Testing:** Should we test S3 direct access from EC2 to measure real-world latency?
2. **Compression Levels:** Default to level 3, or allow configuration?
3. **Chunk Caching:** Should we cache chunks locally or always fetch from S3?
4. **Multi-region:** Support cross-region S3 access?

## Conclusion

v0.3.0 focuses on production readiness through:
- **Performance:** Parallel compression for 3-5x faster conversion
- **Cloud-native:** Direct S3 access eliminates download step
- **Efficiency:** Streaming processing reduces memory usage

These improvements make BAMS3 practical for large-scale genomics workflows in cloud environments, completing the transition from proof-of-concept to production-ready format.

**Target delivery:** February 1, 2026 (2 weeks)
**Status:** Ready to begin implementation

---

**Plan Version:** 1.0
**Date:** 2026-01-15
**Author:** Scott Friedman
