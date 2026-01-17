# CloudNativeFormats Core Library Design

## Overview

A shared library providing common infrastructure for cloud-native scientific data formats. Designed to be reused across BAMS3 (genomic alignments), VCFS3 (variants), RVIDS3 (research video), and future formats (CryoEM, microscopy, etc.).

## Problem Statement

### Current Code Duplication

As we implement multiple cloud-native formats, common patterns emerge:

1. **Storage Abstraction**: Local files vs S3 vs other cloud providers
2. **Chunking Strategy**: Dividing data into queryable chunks
3. **Index Management**: Building and querying indices for selective access
4. **Compression**: Consistent compression with streaming support
5. **Metadata**: Standard format for dataset information
6. **Parallel Processing**: Worker pools and concurrent uploads
7. **Error Handling**: Retries, timeouts, partial failures

**Without a core library**, each format reimplements these, leading to:
- Code duplication (3-5x redundant code)
- Inconsistent behavior across formats
- Difficult to maintain and improve
- Steep learning curve for contributors

**With a core library**, we get:
- Single implementation of common functionality
- Consistent behavior and testing
- Easy to add new formats (focus on domain-specific logic)
- Shared improvements benefit all formats

## Design Goals

1. **Pluggable Architecture**: Easy to swap storage backends, compression algorithms, etc.
2. **Zero-Copy Where Possible**: Minimize memory allocations
3. **Streaming First**: Support streaming reads/writes
4. **Cloud-Native**: Optimized for S3, but not limited to it
5. **Testable**: Mock-friendly interfaces
6. **Performance**: Minimal overhead over direct implementation
7. **Language Support**: Go first, Python bindings later

## Architecture

### Package Structure

```
cloudnative-formats/
├── core/                           # Core library (Go)
│   ├── storage/                    # Storage abstraction
│   │   ├── storage.go              # Interfaces
│   │   ├── local.go                # Local filesystem
│   │   ├── s3.go                   # AWS S3
│   │   └── memory.go               # In-memory (testing)
│   ├── chunk/                      # Chunking strategies
│   │   ├── chunk.go                # Interfaces
│   │   ├── fixed.go                # Fixed-size chunks
│   │   ├── coordinate.go           # Coordinate-based (genomics)
│   │   ├── temporal.go             # Time-based (video)
│   │   └── grid.go                 # 2D/3D grid (imaging)
│   ├── index/                      # Index management
│   │   ├── index.go                # Interfaces
│   │   ├── btree.go                # B-tree index
│   │   ├── interval.go             # Interval tree (genomics)
│   │   └── spatial.go              # Spatial index (imaging)
│   ├── compress/                   # Compression
│   │   ├── compress.go             # Interfaces
│   │   ├── zstd.go                 # Zstandard
│   │   ├── lz4.go                  # LZ4
│   │   └── none.go                 # No compression
│   ├── metadata/                   # Metadata management
│   │   ├── metadata.go             # Standard metadata format
│   │   └── schema.go               # JSON schema validation
│   ├── worker/                     # Parallel processing
│   │   ├── pool.go                 # Worker pool
│   │   └── pipeline.go             # Processing pipeline
│   ├── query/                      # Query engine components
│   │   ├── planner.go              # Query planning
│   │   ├── executor.go             # Query execution
│   │   └── filter.go               # Filtering logic
│   └── util/                       # Utilities
│       ├── buffer.go               # Buffer pools
│       ├── retry.go                # Retry logic
│       └── progress.go             # Progress reporting
├── formats/                        # Format implementations
│   ├── bams3/                      # BAMS3 (genomic alignments)
│   ├── vcfs3/                      # VCFS3 (variants)
│   ├── rvids3/                     # RVIDS3 (research video)
│   └── cryoem/                     # CryoEM (future)
├── bindings/                       # Language bindings
│   └── python/                     # Python bindings (future)
├── docs/                           # Documentation
│   ├── FORMAT_DESIGN_GUIDE.md      # How to design a new format
│   ├── API_REFERENCE.md            # Core API documentation
│   └── EXAMPLES.md                 # Usage examples
└── examples/                       # Example code
    ├── simple_format/              # Minimal format example
    └── benchmarks/                 # Performance benchmarks
```

### Core Interfaces

#### Storage Interface

**Purpose**: Abstract storage backend (local, S3, Azure, GCS, etc.)

```go
package storage

import (
	"context"
	"io"
)

// Storage represents a storage backend
type Storage interface {
	// Open opens a file for reading
	Open(ctx context.Context, path string) (Reader, error)

	// Create creates a file for writing
	Create(ctx context.Context, path string) (Writer, error)

	// List lists files matching a prefix
	List(ctx context.Context, prefix string) ([]FileInfo, error)

	// Delete deletes a file
	Delete(ctx context.Context, path string) error

	// Exists checks if a file exists
	Exists(ctx context.Context, path string) (bool, error)

	// Stat gets file metadata
	Stat(ctx context.Context, path string) (FileInfo, error)
}

// Reader extends io.ReadCloser with range reads
type Reader interface {
	io.ReadCloser
	io.ReaderAt
	io.Seeker

	// ReadRange reads a specific byte range
	ReadRange(offset int64, length int64) ([]byte, error)
}

// Writer extends io.WriteCloser with multipart uploads
type Writer interface {
	io.WriteCloser

	// Abort aborts a partial write (e.g., multipart upload)
	Abort() error
}

// FileInfo provides file metadata
type FileInfo struct {
	Path         string
	Size         int64
	ModTime      time.Time
	ETag         string    // For S3
	StorageClass string    // For S3
}

// Options configures storage behavior
type Options struct {
	// For S3
	Region       string
	Bucket       string
	PartSize     int64  // Multipart upload part size
	Concurrency  int    // Concurrent uploads

	// For local
	RootDir      string

	// Common
	Timeout      time.Duration
	MaxRetries   int
}
```

**Implementations**:

1. **Local Storage** (`local.go`):
   ```go
   type LocalStorage struct {
   	rootDir string
   }

   func NewLocalStorage(rootDir string) *LocalStorage {
   	return &LocalStorage{rootDir: rootDir}
   }
   ```

2. **S3 Storage** (`s3.go`):
   ```go
   type S3Storage struct {
   	client      *s3.Client
   	bucket      string
   	partSize    int64
   	concurrency int
   }

   func NewS3Storage(opts Options) (*S3Storage, error) {
   	// Initialize AWS SDK
   	// Configure multipart uploads
   }
   ```

3. **Memory Storage** (`memory.go`) - for testing:
   ```go
   type MemoryStorage struct {
   	files map[string][]byte
   	mu    sync.RWMutex
   }
   ```

#### Sharding Interface

**Purpose**: Distribute chunks across S3 prefixes for higher aggregate throughput

```go
package storage

import (
	"crypto/sha256"
	"fmt"
)

// Sharding strategy determines how to distribute chunks across prefixes
type ShardingStrategy interface {
	// GetPrefix returns the S3 prefix for a given chunk path
	GetPrefix(chunkPath string) string

	// NumPrefixes returns the total number of prefixes
	NumPrefixes() int

	// IsEnabled returns whether sharding is enabled
	IsEnabled() bool
}

// NoSharding uses sequential prefixes (default, backward compatible)
type NoSharding struct{}

func (n *NoSharding) GetPrefix(chunkPath string) string {
	return ""
}

func (n *NoSharding) NumPrefixes() int {
	return 1
}

func (n *NoSharding) IsEnabled() bool {
	return false
}

// HashSharding distributes chunks using hash-based prefixes
type HashSharding struct {
	Bits int  // 4 bits = 16 prefixes, 6 bits = 64, 8 bits = 256
}

func NewHashSharding(bits int) *HashSharding {
	if bits < 4 || bits > 8 {
		bits = 8  // Default to 8 bits (256 prefixes)
	}
	return &HashSharding{Bits: bits}
}

func (h *HashSharding) GetPrefix(chunkPath string) string {
	// Hash the chunk path
	hash := sha256.Sum256([]byte(chunkPath))

	// Take first N bits
	prefixValue := hash[0] >> (8 - h.Bits)

	// Convert to hex (1 or 2 chars depending on bits)
	if h.Bits <= 4 {
		return fmt.Sprintf("%x", prefixValue)
	}
	return fmt.Sprintf("%02x", prefixValue)
}

func (h *HashSharding) NumPrefixes() int {
	return 1 << h.Bits  // 2^bits
}

func (h *HashSharding) IsEnabled() bool {
	return true
}

// ShardedStorage wraps a storage backend with sharding support
type ShardedStorage struct {
	backend  Storage
	strategy ShardingStrategy
}

func NewShardedStorage(backend Storage, strategy ShardingStrategy) *ShardedStorage {
	return &ShardedStorage{
		backend:  backend,
		strategy: strategy,
	}
}

func (s *ShardedStorage) Open(ctx context.Context, path string) (Reader, error) {
	// Prepend shard prefix
	shardedPath := s.getShardedPath(path)
	return s.backend.Open(ctx, shardedPath)
}

func (s *ShardedStorage) Create(ctx context.Context, path string) (Writer, error) {
	shardedPath := s.getShardedPath(path)
	return s.backend.Create(ctx, shardedPath)
}

func (s *ShardedStorage) getShardedPath(path string) string {
	if !s.strategy.IsEnabled() {
		return path
	}

	prefix := s.strategy.GetPrefix(path)
	if prefix == "" {
		return path
	}

	return prefix + "/" + path
}

// List implementation needs to scan all prefixes
func (s *ShardedStorage) List(ctx context.Context, prefix string) ([]FileInfo, error) {
	if !s.strategy.IsEnabled() {
		return s.backend.List(ctx, prefix)
	}

	// List across all shard prefixes
	var allFiles []FileInfo
	numPrefixes := s.strategy.NumPrefixes()

	for i := 0; i < numPrefixes; i++ {
		shardPrefix := fmt.Sprintf("%02x/%s", i, prefix)
		files, err := s.backend.List(ctx, shardPrefix)
		if err != nil {
			return nil, err
		}
		allFiles = append(allFiles, files...)
	}

	return allFiles, nil
}
```

**Usage in Format Implementation**:
```go
// Without sharding (default)
storage := storage.NewS3Storage(opts)

// With sharding enabled
baseStorage := storage.NewS3Storage(opts)
sharding := storage.NewHashSharding(8)  // 256 prefixes
storage := storage.NewShardedStorage(baseStorage, sharding)

// Both work with same interface
reader, _ := storage.Open(ctx, "chunks/chr1/00000000-01000000.chunk")
```

**When to Enable**:
- Large cohorts (>1,000 samples for VCFS3)
- High query volume (>100 concurrent queries)
- Production pipelines with many users
- Skip for local storage (no benefit)

**Performance**:
- No sharding: 5,500 requests/sec (single prefix limit)
- 8-bit sharding: 1.4M requests/sec (256 prefixes × 5,500)
- 100× throughput increase for concurrent workloads

#### Chunking Interface

**Purpose**: Define how data is divided into chunks

```go
package chunk

import (
	"context"
	"io"
)

// Strategy defines how to chunk data
type Strategy interface {
	// Plan determines chunk boundaries for input data
	Plan(ctx context.Context, input io.Reader, metadata map[string]interface{}) ([]ChunkInfo, error)

	// CreateWriter creates a writer for a specific chunk
	CreateWriter(ctx context.Context, chunkID string) (io.WriteCloser, error)

	// CreateReader creates a reader for a specific chunk
	CreateReader(ctx context.Context, chunkID string) (io.ReadCloser, error)
}

// ChunkInfo describes a chunk
type ChunkInfo struct {
	ID       string                 // Unique chunk identifier
	Size     int64                  // Size in bytes
	Metadata map[string]interface{} // Format-specific metadata
}

// Coordinator manages chunking for a dataset
type Coordinator struct {
	strategy Strategy
	storage  storage.Storage
	compress compress.Compressor
}

func NewCoordinator(strategy Strategy, storage storage.Storage, compress compress.Compressor) *Coordinator {
	return &Coordinator{
		strategy: strategy,
		storage:  storage,
		compress: compress,
	}
}
```

**Implementations**:

1. **Fixed-Size Chunks** (`fixed.go`):
   ```go
   type FixedSizeStrategy struct {
   	chunkSize int64
   }
   ```

2. **Coordinate-Based Chunks** (`coordinate.go`) - for genomics:
   ```go
   type CoordinateStrategy struct {
   	chunkSize    int64  // e.g., 1Mbp
   	sortByCoord  bool
   }
   ```

3. **Temporal Chunks** (`temporal.go`) - for video:
   ```go
   type TemporalStrategy struct {
   	duration time.Duration  // e.g., 30 seconds
   }
   ```

4. **Grid Chunks** (`grid.go`) - for imaging:
   ```go
   type GridStrategy struct {
   	tileWidth  int
   	tileHeight int
   	tileDepth  int  // For 3D volumes
   }
   ```

#### Index Interface

**Purpose**: Fast lookups for selective access

```go
package index

import (
	"context"
)

// Index provides fast lookups
type Index interface {
	// Build creates an index from chunks
	Build(ctx context.Context, chunks []chunk.ChunkInfo) error

	// Query finds chunks matching a query
	Query(ctx context.Context, query Query) ([]chunk.ChunkInfo, error)

	// Save serializes the index
	Save(ctx context.Context, writer io.Writer) error

	// Load deserializes the index
	Load(ctx context.Context, reader io.Reader) error
}

// Query represents a query against an index
type Query interface {
	// Format-specific query implementations
}

// Builder constructs an index incrementally
type Builder interface {
	Add(key interface{}, value chunk.ChunkInfo) error
	Build() (Index, error)
}
```

**Implementations**:

1. **B-Tree Index** (`btree.go`) - general purpose:
   ```go
   type BTreeIndex struct {
   	tree *btree.BTree
   }
   ```

2. **Interval Tree** (`interval.go`) - for genomic coordinates:
   ```go
   type IntervalIndex struct {
   	intervals map[string]*IntervalTree  // Per chromosome
   }

   type IntervalQuery struct {
   	Chromosome string
   	Start      int64
   	End        int64
   }
   ```

3. **Spatial Index** (`spatial.go`) - for imaging:
   ```go
   type SpatialIndex struct {
   	rtree *rtree.RTree
   }

   type SpatialQuery struct {
   	X, Y, Z      int
   	Width, Height, Depth int
   }
   ```

#### Compression Interface

**Purpose**: Pluggable compression algorithms

```go
package compress

import (
	"io"
)

// Compressor provides compression/decompression
type Compressor interface {
	// Compress compresses data
	Compress(src io.Reader, dst io.Writer) error

	// Decompress decompresses data
	Decompress(src io.Reader, dst io.Writer) error

	// Extension returns file extension (e.g., ".zst")
	Extension() string

	// Level returns compression level
	Level() int
}

// Options configures compression
type Options struct {
	Algorithm string  // "zstd", "lz4", "none"
	Level     int     // Compression level (1-22 for zstd)
}

func New(opts Options) (Compressor, error) {
	switch opts.Algorithm {
	case "zstd":
		return NewZstd(opts.Level), nil
	case "lz4":
		return NewLZ4(opts.Level), nil
	case "none":
		return NewNone(), nil
	default:
		return nil, fmt.Errorf("unknown algorithm: %s", opts.Algorithm)
	}
}
```

**Implementations**:

1. **Zstandard** (`zstd.go`):
   ```go
   type ZstdCompressor struct {
   	level int
   }
   ```

2. **LZ4** (`lz4.go`):
   ```go
   type LZ4Compressor struct {
   	level int
   }
   ```

3. **None** (`none.go`) - passthrough:
   ```go
   type NoneCompressor struct{}
   ```

#### Metadata Interface

**Purpose**: Standard dataset metadata format

```go
package metadata

import (
	"time"
)

// Metadata describes a dataset
type Metadata struct {
	Format      string                 `json:"format"`       // e.g., "BAMS3", "VCFS3"
	Version     string                 `json:"version"`      // Format version
	Created     time.Time              `json:"created"`
	DatasetID   string                 `json:"dataset_id"`
	Description string                 `json:"description"`

	// Storage information
	Storage     StorageInfo            `json:"storage"`

	// Chunking information
	Chunks      ChunkingInfo           `json:"chunks"`

	// Compression information
	Compression CompressionInfo        `json:"compression"`

	// Format-specific metadata
	Custom      map[string]interface{} `json:"custom"`
}

type StorageInfo struct {
	Backend   string `json:"backend"`    // "local", "s3", "azure", "gcs"
	Location  string `json:"location"`   // Path or URI
	TotalSize int64  `json:"total_size"` // Total bytes
}

type ChunkingInfo struct {
	Strategy   string `json:"strategy"`    // "fixed", "coordinate", "temporal", etc.
	NumChunks  int    `json:"num_chunks"`
	ChunkSize  int64  `json:"chunk_size"`  // Average or fixed size
}

type CompressionInfo struct {
	Algorithm string  `json:"algorithm"`
	Level     int     `json:"level"`
	Ratio     float64 `json:"ratio"`  // Compression ratio
}

// Manager handles metadata persistence
type Manager struct {
	storage storage.Storage
}

func (m *Manager) Save(ctx context.Context, path string, meta Metadata) error {
	// Serialize to JSON
	// Write to storage
}

func (m *Manager) Load(ctx context.Context, path string) (Metadata, error) {
	// Read from storage
	// Deserialize from JSON
}
```

#### Worker Pool Interface

**Purpose**: Parallel processing with backpressure

```go
package worker

import (
	"context"
	"sync"
)

// Pool manages a pool of workers
type Pool struct {
	workers   int
	taskQueue chan Task
	wg        sync.WaitGroup
}

// Task represents a unit of work
type Task interface {
	Execute(ctx context.Context) error
}

func NewPool(workers int) *Pool {
	return &Pool{
		workers:   workers,
		taskQueue: make(chan Task, workers*2),
	}
}

func (p *Pool) Submit(task Task) {
	p.taskQueue <- task
}

func (p *Pool) Wait() error {
	close(p.taskQueue)
	p.wg.Wait()
	return nil
}

// Pipeline processes data through multiple stages
type Pipeline struct {
	stages []Stage
}

type Stage interface {
	Process(ctx context.Context, input <-chan interface{}, output chan<- interface{}) error
}
```

### Query Engine Components

#### Query Planner

**Purpose**: Optimize query execution

```go
package query

import (
	"context"
)

// Planner creates execution plans
type Planner struct {
	index    index.Index
	metadata metadata.Metadata
}

// Plan creates an execution plan for a query
type Plan struct {
	Chunks      []chunk.ChunkInfo  // Chunks to access
	Filters     []Filter           // Filters to apply
	Order       []OrderBy          // Ordering
	Limit       int                // Result limit
	EstimatedCost Cost             // Cost estimate
}

type Cost struct {
	BytesToDownload int64
	ChunksToAccess  int
	EstimatedTime   time.Duration
}

func (p *Planner) Plan(ctx context.Context, query interface{}) (*Plan, error) {
	// Use index to identify relevant chunks
	// Estimate cost
	// Optimize execution order
}
```

#### Query Executor

**Purpose**: Execute query plans efficiently

```go
package query

import (
	"context"
)

// Executor executes query plans
type Executor struct {
	storage  storage.Storage
	pool     *worker.Pool
	compress compress.Compressor
}

func (e *Executor) Execute(ctx context.Context, plan *Plan) (Results, error) {
	// Download chunks in parallel
	// Apply filters
	// Aggregate results
}

type Results interface {
	Next() bool
	Value() interface{}
	Err() error
	Close() error
}
```

## Format Implementation Pattern

### Minimal Format Example

```go
package simpleformat

import (
	"github.com/cloudnative-formats/core/storage"
	"github.com/cloudnative-formats/core/chunk"
	"github.com/cloudnative-formats/core/compress"
	"github.com/cloudnative-formats/core/metadata"
)

// SimpleFormat is a minimal format implementation
type SimpleFormat struct {
	storage   storage.Storage
	chunking  chunk.Strategy
	compress  compress.Compressor
	meta      metadata.Manager
}

func New(storageURI string) (*SimpleFormat, error) {
	// Parse URI to determine storage backend
	stor, err := storage.FromURI(storageURI)
	if err != nil {
		return nil, err
	}

	// Use fixed-size chunking
	chunking := chunk.NewFixedSize(10 * 1024 * 1024) // 10MB chunks

	// Use zstd compression
	compress, _ := compress.New(compress.Options{
		Algorithm: "zstd",
		Level:     3,
	})

	return &SimpleFormat{
		storage:  stor,
		chunking: chunking,
		compress: compress,
		meta:     metadata.NewManager(stor),
	}, nil
}

func (f *SimpleFormat) Write(ctx context.Context, data io.Reader) error {
	// Plan chunks
	chunks, err := f.chunking.Plan(ctx, data, nil)
	if err != nil {
		return err
	}

	// Write chunks in parallel
	pool := worker.NewPool(8)
	for _, chunkInfo := range chunks {
		pool.Submit(&writeTask{
			format: f,
			chunk:  chunkInfo,
		})
	}
	pool.Wait()

	// Save metadata
	meta := metadata.Metadata{
		Format:  "SimpleFormat",
		Version: "0.1.0",
		// ... fill in metadata
	}
	return f.meta.Save(ctx, "_metadata.json", meta)
}

func (f *SimpleFormat) Read(ctx context.Context, output io.Writer) error {
	// Load metadata
	meta, err := f.meta.Load(ctx, "_metadata.json")
	if err != nil {
		return err
	}

	// Read chunks
	for i := 0; i < meta.Chunks.NumChunks; i++ {
		reader, err := f.chunking.CreateReader(ctx, fmt.Sprintf("chunk_%06d", i))
		if err != nil {
			return err
		}
		defer reader.Close()

		// Decompress and write
		if err := f.compress.Decompress(reader, output); err != nil {
			return err
		}
	}

	return nil
}
```

## Migration Plan

### Phase 1: Extract Core Library (1 week)

1. **Create repository structure**:
   ```bash
   mkdir -p cloudnative-formats/core/{storage,chunk,index,compress,metadata,worker,query,util}
   ```

2. **Extract storage abstraction from BAMS3**:
   - Move `storage.go` to core library
   - Keep S3 and local implementations
   - Update BAMS3 to import from core

3. **Extract compression**:
   - Move `compress.go` to core library
   - Generic zstd wrapper

4. **Basic tests**:
   - Unit tests for each component
   - Integration tests with local storage

### Phase 2: Migrate BAMS3 (1 week)

1. **Update imports**:
   ```go
   import (
   	"github.com/cloudnative-formats/core/storage"
   	"github.com/cloudnative-formats/core/compress"
   )
   ```

2. **Use core chunking**:
   - Implement `CoordinateStrategy` in core
   - Migrate BAMS3 to use it

3. **Use core metadata**:
   - Define BAMS3-specific metadata in `custom` field
   - Use standard metadata manager

4. **Validate**:
   - Run existing BAMS3 tests
   - Verify no performance regression

### Phase 3: Implement VCFS3 with Core (2 weeks)

1. **Use core from start**:
   - Storage: Use core `storage.Storage`
   - Chunking: Implement 2D chunking strategy in core
   - Index: Use `IntervalIndex` from core
   - Compression: Use core `compress.Compressor`

2. **Add VCFS3-specific logic**:
   - Sparse genotype encoding (format-specific)
   - VCF parsing and export (format-specific)
   - Sample manifest (format-specific)

3. **Demonstrate reuse**:
   - Document how much code is shared vs custom
   - Show reduced development time

### Phase 4: Documentation and Examples (1 week)

1. **API Documentation**:
   - godoc for all core packages
   - Interface documentation with examples

2. **Format Design Guide** (`docs/FORMAT_DESIGN_GUIDE.md`):
   - Step-by-step guide to creating a new format
   - Decision tree for choosing chunking strategies
   - Performance considerations

3. **Examples**:
   - Simple format (minimal example)
   - Complex format (using all features)
   - Benchmarks comparing approaches

### Phase 5: Python Bindings (2 weeks, optional)

1. **CGo bindings** or **gRPC server**:
   - Expose core library to Python
   - Python-friendly API

2. **Python package**:
   ```python
   from cloudnative_formats import Storage, Chunking, Compression

   storage = Storage.from_uri("s3://bucket/dataset")
   format = MyFormat(storage)
   format.write(data)
   ```

## Benefits

### Code Reuse

**Before Core Library** (BAMS3 + VCFS3):
- Storage code: 500 lines × 2 = 1,000 lines
- Compression: 200 lines × 2 = 400 lines
- Worker pool: 150 lines × 2 = 300 lines
- **Total**: 1,700 lines duplicated

**After Core Library**:
- Core library: 2,000 lines (one implementation)
- BAMS3 specific: 500 lines
- VCFS3 specific: 600 lines
- **Total**: 3,100 lines (vs 4,700 lines)
- **Savings**: 34% reduction

### Consistency

All formats get:
- Same retry logic
- Same error handling
- Same progress reporting
- Same performance optimizations
- Same bug fixes

### Development Speed

**Time to implement new format**:
- **Without core**: 2-4 weeks (reimplementing everything)
- **With core**: 3-5 days (focus on domain logic)

**Example**: RVIDS3 implementation
- Video segmentation: 2 days (domain-specific)
- Annotation indexing: 1 day (domain-specific)
- Storage/compression/chunking: 0 days (core library)
- **Total**: 3 days vs 2-3 weeks

### Testing

**Core library tests** (once):
- Storage backends: 200 test cases
- Compression: 50 test cases
- Chunking: 100 test cases
- **Total**: 350 test cases

**Format tests** (per format):
- Domain-specific logic: 100 test cases
- Integration with core: 20 test cases
- **Total per format**: 120 test cases

**Result**: Higher test coverage with less effort

## Performance Considerations

### Zero-Copy Operations

Use `io.Reader`/`io.Writer` interfaces to avoid copying:

```go
// Good: zero-copy streaming
func (s *S3Storage) Read(ctx context.Context, path string) (io.ReadCloser, error) {
	// Returns S3 stream directly
}

// Bad: copies entire file into memory
func (s *S3Storage) ReadAll(ctx context.Context, path string) ([]byte, error) {
	// Avoid this pattern
}
```

### Buffer Pooling

Reuse buffers to reduce GC pressure:

```go
package util

import "sync"

var bufferPool = sync.Pool{
	New: func() interface{} {
		return make([]byte, 64*1024)
	},
}

func GetBuffer() []byte {
	return bufferPool.Get().([]byte)
}

func PutBuffer(buf []byte) {
	bufferPool.Put(buf)
}
```

### Parallel Processing

Use worker pools for I/O-bound operations:

```go
// Good: parallel chunk downloads
pool := worker.NewPool(10)
for _, chunk := range chunks {
	pool.Submit(&downloadTask{chunk: chunk})
}
pool.Wait()

// Bad: sequential downloads
for _, chunk := range chunks {
	downloadChunk(chunk)
}
```

## Testing Strategy

### Unit Tests

Each core component has comprehensive unit tests:

```go
package storage

func TestS3Storage_ReadRange(t *testing.T) {
	// Test S3 range requests
}

func TestLocalStorage_Concurrent(t *testing.T) {
	// Test concurrent access
}
```

### Integration Tests

Test interactions between components:

```go
package integration

func TestChunkingWithCompression(t *testing.T) {
	storage := storage.NewMemory()
	compress := compress.NewZstd(3)
	chunking := chunk.NewFixed(1024)

	// Test full pipeline
}
```

### Format Tests

Each format tests core integration:

```go
package bams3

func TestBAMS3_WithCoreLibrary(t *testing.T) {
	// Verify BAMS3 works with core components
}
```

### Performance Tests

Benchmark core operations:

```go
func BenchmarkS3_RangeRequest(b *testing.B) {
	// Measure S3 range request performance
}

func BenchmarkZstd_Compress(b *testing.B) {
	// Measure compression performance
}
```

## Documentation

### Format Design Guide

`docs/FORMAT_DESIGN_GUIDE.md` covers:

1. **Choosing a chunking strategy**:
   - Fixed-size: Simple, general-purpose
   - Coordinate-based: Genomic data
   - Temporal: Time-series data
   - Grid: Multi-dimensional imaging

2. **Choosing an index type**:
   - B-tree: Sorted keys
   - Interval tree: Range queries
   - Spatial: Multi-dimensional queries

3. **Compression tradeoffs**:
   - zstd: Best compression, medium speed
   - lz4: Fast, lower compression
   - none: No overhead, larger files

4. **Metadata design**:
   - Standard fields (format, version, created)
   - Custom fields (format-specific)

5. **Performance optimization**:
   - Chunk size selection
   - Parallel processing
   - Index design

### API Reference

`docs/API_REFERENCE.md` provides:

- Interface documentation
- Usage examples
- Best practices
- Common patterns

## Example: CryoEM Format

Using the core library to implement a CryoEM format:

```go
package cryoem

import (
	"github.com/cloudnative-formats/core/storage"
	"github.com/cloudnative-formats/core/chunk"
	"github.com/cloudnative-formats/core/compress"
	"github.com/cloudnative-formats/core/index"
)

type CryoEMFormat struct {
	storage   storage.Storage
	chunking  *chunk.GridStrategy    // 3D tiles
	compress  compress.Compressor     // LZ4 for speed
	index     *index.SpatialIndex    // 3D spatial index
}

func New(uri string) (*CryoEMFormat, error) {
	stor, _ := storage.FromURI(uri)

	return &CryoEMFormat{
		storage: stor,
		chunking: chunk.NewGrid(
			256,  // tile width
			256,  // tile height
			32,   // tile depth (z-stack)
		),
		compress: compress.NewLZ4(1),  // Fast decompression
		index:    index.NewSpatial(),
	}, nil
}

func (f *CryoEMFormat) Write(ctx context.Context, volume []byte, dims [3]int) error {
	// Domain-specific: 3D volume slicing
	tiles := f.sliceVolume(volume, dims)

	// Core library: parallel write with compression
	pool := worker.NewPool(16)
	for _, tile := range tiles {
		pool.Submit(&writeTile{
			storage:  f.storage,
			compress: f.compress,
			tile:     tile,
		})
	}
	return pool.Wait()
}

func (f *CryoEMFormat) Query(ctx context.Context, x, y, z, width, height, depth int) ([]byte, error) {
	// Core library: spatial index query
	query := &index.SpatialQuery{
		X: x, Y: y, Z: z,
		Width: width, Height: height, Depth: depth,
	}

	chunks, err := f.index.Query(ctx, query)
	if err != nil {
		return nil, err
	}

	// Domain-specific: reconstruct volume from tiles
	return f.reconstructVolume(ctx, chunks)
}
```

**Result**: CryoEM format implemented in ~500 lines vs ~2,000 lines without core library

## Summary

The core library provides:

- **Storage abstraction**: Local, S3, Azure, GCS with consistent interface
- **Chunking strategies**: Fixed, coordinate, temporal, grid
- **Index types**: B-tree, interval tree, spatial index
- **Compression**: zstd, LZ4, pluggable
- **Parallel processing**: Worker pools, pipelines
- **Metadata**: Standard format with custom extensions
- **Query engine**: Planning and execution

**Benefits**:
- 34% code reduction
- 5x faster format development
- Consistent behavior across formats
- Shared improvements and bug fixes
- Lower learning curve

**Status**: Design complete, ready for implementation

**Next Steps**:
1. Extract core library from BAMS3 (1 week)
2. Migrate BAMS3 to use core (1 week)
3. Implement VCFS3 with core (2 weeks)
4. Document format design guide (1 week)
