# S3 direct integration for cloud-native workflows

**Labels:** enhancement, cloud, priority-high
**Priority:** High

## Goal
Enable direct read/write from S3 without local storage, demonstrating true cloud-native architecture.

## Implementation Plan

### 1. S3 Writer Integration
**Path:** `pkg/bams3/s3_writer.go` (new file)

```go
type S3Writer struct {
    client *s3.Client
    bucket string
    prefix string
    uploader *manager.Uploader
}
```

**Features:**
- [ ] Multipart upload for large chunks
- [ ] Parallel chunk uploads (leverage existing workers)
- [ ] Progress tracking per chunk
- [ ] Retry logic with exponential backoff
- [ ] Metadata and header upload

### 2. S3 Reader Integration
**Path:** `pkg/bams3/s3_reader.go` (new file)

```go
type S3Reader struct {
    client *s3.Client
    bucket string
    prefix string
}
```

**Features:**
- [ ] Range requests for selective chunk downloads
- [ ] Parallel chunk downloads
- [ ] Metadata/header caching
- [ ] Index caching for fast queries

### 3. CLI Integration

**Convert to S3:**
```bash
bwa mem ref.fa reads.fq | bams3 convert --stdin s3://bucket/sample.bams3
```

**Query from S3:**
```bash
bams3 query s3://bucket/sample.bams3 chr1:1000000-2000000
bams3 stats s3://bucket/sample.bams3
```

### 4. Performance Optimizations
- [ ] Chunk upload concurrency (default: 2x workers)
- [ ] Part size tuning (default: 10MB)
- [ ] Connection pooling
- [ ] Compression before upload
- [ ] Same-region optimization (free transfer)

## Configuration

Add S3 flags to convert command:
```go
--s3-part-size     // Default: 10MB
--s3-concurrency   // Default: workers * 2
--s3-region        // Auto-detect if possible
```

## Testing

### Unit Tests
- [ ] Multipart upload logic
- [ ] Range request handling
- [ ] Error handling and retries

### Integration Tests
- [ ] Upload small dataset to S3
- [ ] Query from S3
- [ ] Large file upload (>5GB)
- [ ] Cross-region transfer

### Performance Tests
- [ ] Compare upload speed vs local disk
- [ ] Measure query latency from S3
- [ ] Test with different part sizes

## Acceptance Criteria
- [ ] Can stream directly to S3 without local storage
- [ ] Can query directly from S3
- [ ] Performance comparable to local disk (same region)
- [ ] Proper error handling and retries
- [ ] Documentation with examples
- [ ] Cost estimates in README

## Dependencies
- AWS SDK for Go v2: `github.com/aws/aws-sdk-go-v2/service/s3`
- S3 Transfer Manager: `github.com/aws/aws-sdk-go-v2/feature/s3/manager`

## Notes
This is the core value proposition of cloud-native format.
