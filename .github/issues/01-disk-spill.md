# Complete disk spill implementation (external sort)

**Labels:** enhancement, performance
**Priority:** Medium

## Goal
Enable processing of datasets that exceed available RAM by implementing external sort with disk spill.

## Current Status
Infrastructure is in place:
- ✅ Temporary directory creation (`/tmp/bams3-spill-*`)
- ✅ SpillFile tracking structure
- ✅ Automatic cleanup on success/failure
- ✅ Configuration flag (`SpillToDisk`, default: true)

## Remaining Work

### 1. Implement Spill Logic
- [ ] Add memory pressure detection in `flushSortBuffer`
- [ ] Write sorted buffer to temporary file when spilling
- [ ] Use msgpack or similar for efficient serialization
- [ ] Track spill files in `spillFiles` array
- [ ] Update `SpillsToDisK` counter in stats

### 2. Implement K-Way Merge
- [ ] Create `SpillReader` for reading from spill files
- [ ] Implement min-heap for k-way merge
- [ ] Stream merged results to chunk writer
- [ ] Handle reference boundaries correctly

### 3. Update Progress Reporting
- [ ] Show spill count in progress output
- [ ] Display spill count in final stats

### 4. Testing
- [ ] Test with dataset requiring 2 spills
- [ ] Test with dataset requiring 10+ spills
- [ ] Verify correctness (sorted output, all reads present)
- [ ] Benchmark performance vs in-memory sort

## Acceptance Criteria
- [ ] Can process 100GB+ input with 8GB sort buffer
- [ ] Spilled data correctly merged in coordinate order
- [ ] Performance within 2x of in-memory sort
- [ ] Temporary files cleaned up on success and failure
- [ ] Progress reporting shows spill activity

## Related Code
- `pkg/bams3/stream_converter.go:131-152` (SpillFile struct)
- `pkg/bams3/stream_converter.go:392-470` (flushSortBuffer)
- `pkg/bams3/stream_converter.go:761-774` (cleanup)

## Notes
Current incremental flushing handles most real-world cases. This is needed for extreme datasets (whole genomes at high coverage).
