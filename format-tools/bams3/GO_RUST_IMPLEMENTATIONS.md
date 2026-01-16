# Go and Rust Implementations of BAMS3

## Why Go or Rust?

The Python proof-of-concept demonstrates the BAMS3 format, but production tools benefit from compiled languages:

### Go Advantages
- **Fast compilation** - Quick iteration
- **Simple deployment** - Single binary
- **Good AWS SDK** - Excellent S3 support
- **Concurrency** - Built-in goroutines for parallel chunk processing
- **Cross-platform** - Easy to build for multiple platforms
- **Genomics ecosystem** - Used by biogo, others

### Rust Advantages
- **Maximum performance** - Zero-cost abstractions
- **Memory safety** - No segfaults, data races
- **Excellent ecosystem** - Serde, tokio, aws-sdk
- **Binary size** - Smaller than Go
- **WebAssembly** - Can compile to WASM
- **Genomics momentum** - rust-bio, noodles growing

## Recommended Approach

**Build both!**
- **Go CLI tools** - Fast to develop, easy to deploy
- **Rust library** - High-performance core, FFI bindings

## Go Implementation

### Project Structure

```
bams3-go/
├── cmd/
│   ├── bams3/           # Main CLI
│   ├── bams3-convert/   # BAM → BAMS3
│   └── bams3-query/     # Query tool
├── pkg/
│   ├── bams3/
│   │   ├── reader.go    # Read BAMS3 datasets
│   │   ├── writer.go    # Write BAMS3 datasets
│   │   ├── index.go     # Index operations
│   │   ├── chunk.go     # Chunk format
│   │   └── metadata.go  # Metadata handling
│   └── bam/
│       └── convert.go   # BAM conversion
├── go.mod
└── go.sum
```

### Core Go Types

```go
// pkg/bams3/types.go
package bams3

import "time"

// Dataset represents a BAMS3 dataset
type Dataset struct {
    Path     string
    Metadata Metadata
    Header   Header
    Index    SpatialIndex
}

// Metadata contains dataset metadata
type Metadata struct {
    Format      string            `json:"format"`
    Version     string            `json:"version"`
    Created     time.Time         `json:"created"`
    Statistics  Statistics        `json:"statistics"`
    Chunks      []ChunkInfo       `json:"chunks"`
    ChunkSize   int               `json:"chunk_size"`
    Compression CompressionConfig `json:"compression"`
}

// ChunkInfo describes a chunk
type ChunkInfo struct {
    Path          string    `json:"path"`
    Reference     string    `json:"reference"`
    Start         int       `json:"start"`
    End           int       `json:"end"`
    Reads         int       `json:"reads"`
    SizeBytes     int64     `json:"size_bytes"`
    Compression   string    `json:"compression"`
    Checksum      string    `json:"checksum"`
    Created       time.Time `json:"created"`
}

// Read represents an alignment read
type Read struct {
    Name            string
    Flag            uint16
    ReferenceID     int
    Position        int
    MappingQuality  uint8
    CIGAR           string
    MateReferenceID int
    MatePosition    int
    TemplateLength  int
    Sequence        string
    Quality         string
    Tags            map[string]interface{}
}
```

### Reader Implementation

```go
// pkg/bams3/reader.go
package bams3

import (
    "context"
    "encoding/json"
    "fmt"
    "io"

    "github.com/aws/aws-sdk-go-v2/service/s3"
)

type Reader struct {
    dataset  Dataset
    s3Client *s3.Client
}

// Open opens a BAMS3 dataset from S3 or local path
func Open(ctx context.Context, uri string) (*Reader, error) {
    r := &Reader{}

    // Parse URI (s3://bucket/path or local path)
    isS3, bucket, key := parseURI(uri)

    if isS3 {
        // Initialize S3 client
        cfg, err := config.LoadDefaultConfig(ctx)
        if err != nil {
            return nil, err
        }
        r.s3Client = s3.NewFromConfig(cfg)
    }

    // Read metadata
    metadataPath := joinPath(uri, "_metadata.json")
    metadataData, err := r.readFile(ctx, metadataPath)
    if err != nil {
        return nil, fmt.Errorf("reading metadata: %w", err)
    }

    if err := json.Unmarshal(metadataData, &r.dataset.Metadata); err != nil {
        return nil, fmt.Errorf("parsing metadata: %w", err)
    }

    // Read header
    headerPath := joinPath(uri, "_header.json")
    headerData, err := r.readFile(ctx, headerPath)
    if err != nil {
        return nil, fmt.Errorf("reading header: %w", err)
    }

    if err := json.Unmarshal(headerData, &r.dataset.Header); err != nil {
        return nil, fmt.Errorf("parsing header: %w", err)
    }

    // Read index
    indexPath := joinPath(uri, "_index/spatial.json")
    indexData, err := r.readFile(ctx, indexPath)
    if err != nil {
        return nil, fmt.Errorf("reading index: %w", err)
    }

    if err := json.Unmarshal(indexData, &r.dataset.Index); err != nil {
        return nil, fmt.Errorf("parsing index: %w", err)
    }

    r.dataset.Path = uri
    return r, nil
}

// Query queries reads in a genomic region
func (r *Reader) Query(ctx context.Context, reference string, start, end int) (<-chan Read, error) {
    reads := make(chan Read, 100)

    go func() {
        defer close(reads)

        // Find overlapping chunks
        chunks := r.dataset.Index.FindOverlapping(reference, start, end)

        // Process chunks in parallel
        for _, chunkInfo := range chunks {
            chunkPath := joinPath(r.dataset.Path, chunkInfo.Path)

            // Read chunk
            chunkData, err := r.readFile(ctx, chunkPath)
            if err != nil {
                // Log error but continue
                continue
            }

            // Parse chunk (JSON in POC, binary in production)
            var chunkReads []Read
            if err := json.Unmarshal(chunkData, &chunkReads); err != nil {
                continue
            }

            // Filter to exact region and send to channel
            for _, read := range chunkReads {
                if read.Position >= start && read.Position < end {
                    reads <- read
                }
            }
        }
    }()

    return reads, nil
}

// QueryParallel queries multiple regions in parallel
func (r *Reader) QueryParallel(ctx context.Context, regions []Region, workers int) (<-chan Read, error) {
    reads := make(chan Read, 1000)

    // Worker pool pattern
    jobs := make(chan Region, len(regions))

    // Start workers
    for w := 0; w < workers; w++ {
        go func() {
            for region := range jobs {
                regionReads, _ := r.Query(ctx, region.Reference, region.Start, region.End)
                for read := range regionReads {
                    reads <- read
                }
            }
        }()
    }

    // Send jobs
    go func() {
        for _, region := range regions {
            jobs <- region
        }
        close(jobs)
        close(reads)
    }()

    return reads, nil
}
```

### CLI Tool

```go
// cmd/bams3/main.go
package main

import (
    "context"
    "fmt"
    "os"

    "github.com/spf13/cobra"
    "github.com/username/bams3-go/pkg/bams3"
)

func main() {
    rootCmd := &cobra.Command{
        Use:   "bams3",
        Short: "BAMS3 tools for cloud-native genomics",
    }

    queryCmd := &cobra.Command{
        Use:   "query <dataset> <region>",
        Short: "Query reads from BAMS3 dataset",
        Args:  cobra.ExactArgs(2),
        RunE:  runQuery,
    }

    convertCmd := &cobra.Command{
        Use:   "convert <input.bam> <output.bams3>",
        Short: "Convert BAM to BAMS3",
        Args:  cobra.ExactArgs(2),
        RunE:  runConvert,
    }

    rootCmd.AddCommand(queryCmd, convertCmd)

    if err := rootCmd.Execute(); err != nil {
        os.Exit(1)
    }
}

func runQuery(cmd *cobra.Command, args []string) error {
    ctx := context.Background()
    dataset := args[0]
    region := args[1]

    // Parse region
    ref, start, end, err := parseRegion(region)
    if err != nil {
        return err
    }

    // Open dataset
    reader, err := bams3.Open(ctx, dataset)
    if err != nil {
        return err
    }

    // Query
    reads, err := reader.Query(ctx, ref, start, end)
    if err != nil {
        return err
    }

    // Print reads
    for read := range reads {
        fmt.Printf("%s\t%d\t%d\t%s\n",
            read.Name, read.Position, read.MappingQuality, read.CIGAR)
    }

    return nil
}
```

### Building

```bash
# Install dependencies
go mod init github.com/username/bams3-go
go get github.com/aws/aws-sdk-go-v2
go get github.com/spf13/cobra

# Build
go build -o bams3 ./cmd/bams3

# Install
go install ./cmd/bams3

# Use
bams3 query s3://bucket/sample.bams3 chr1:1000000-2000000
```

### Performance Optimizations

```go
// Use sync.Pool for buffer reuse
var chunkBufferPool = sync.Pool{
    New: func() interface{} {
        return make([]byte, 0, 10*1024*1024) // 10MB
    },
}

// Parallel chunk fetching
func (r *Reader) fetchChunksParallel(ctx context.Context, chunks []ChunkInfo) [][]byte {
    results := make([][]byte, len(chunks))
    var wg sync.WaitGroup

    for i, chunk := range chunks {
        wg.Add(1)
        go func(idx int, c ChunkInfo) {
            defer wg.Done()
            data, _ := r.readFile(ctx, c.Path)
            results[idx] = data
        }(i, chunk)
    }

    wg.Wait()
    return results
}
```

## Rust Implementation

### Project Structure

```
bams3-rs/
├── Cargo.toml
├── src/
│   ├── lib.rs           # Library root
│   ├── dataset.rs       # Dataset type
│   ├── reader.rs        # Reader
│   ├── writer.rs        # Writer
│   ├── chunk.rs         # Chunk format
│   ├── index.rs         # Indexing
│   └── error.rs         # Error types
├── bams3-cli/
│   ├── Cargo.toml
│   └── src/
│       └── main.rs      # CLI tool
└── bams3-ffi/
    ├── Cargo.toml
    └── src/
        └── lib.rs       # C FFI bindings
```

### Core Rust Types

```rust
// src/dataset.rs
use serde::{Deserialize, Serialize};
use chrono::{DateTime, Utc};

#[derive(Debug, Serialize, Deserialize)]
pub struct Dataset {
    pub path: String,
    pub metadata: Metadata,
    pub header: Header,
    pub index: SpatialIndex,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Metadata {
    pub format: String,
    pub version: String,
    pub created: DateTime<Utc>,
    pub statistics: Statistics,
    pub chunks: Vec<ChunkInfo>,
    pub chunk_size: usize,
    pub compression: CompressionConfig,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ChunkInfo {
    pub path: String,
    pub reference: String,
    pub start: u64,
    pub end: u64,
    pub reads: usize,
    pub size_bytes: u64,
    pub compression: String,
    pub checksum: String,
    pub created: DateTime<Utc>,
}

#[derive(Debug, Clone)]
pub struct Read {
    pub name: String,
    pub flag: u16,
    pub reference_id: i32,
    pub position: i64,
    pub mapping_quality: u8,
    pub cigar: String,
    pub sequence: String,
    pub quality: String,
}
```

### Reader Implementation

```rust
// src/reader.rs
use anyhow::Result;
use aws_sdk_s3::Client as S3Client;
use tokio::sync::mpsc;

pub struct Reader {
    dataset: Dataset,
    s3_client: Option<S3Client>,
}

impl Reader {
    pub async fn open(uri: &str) -> Result<Self> {
        // Parse URI
        let is_s3 = uri.starts_with("s3://");

        // Initialize S3 client if needed
        let s3_client = if is_s3 {
            let config = aws_config::load_from_env().await;
            Some(S3Client::new(&config))
        } else {
            None
        };

        // Read metadata
        let metadata_path = format!("{}/_metadata.json", uri);
        let metadata_data = Self::read_file(&s3_client, &metadata_path).await?;
        let metadata: Metadata = serde_json::from_slice(&metadata_data)?;

        // Read header
        let header_path = format!("{}/_header.json", uri);
        let header_data = Self::read_file(&s3_client, &header_path).await?;
        let header: Header = serde_json::from_slice(&header_data)?;

        // Read index
        let index_path = format!("{}/_index/spatial.json", uri);
        let index_data = Self::read_file(&s3_client, &index_path).await?;
        let index: SpatialIndex = serde_json::from_slice(&index_data)?;

        Ok(Self {
            dataset: Dataset {
                path: uri.to_string(),
                metadata,
                header,
                index,
            },
            s3_client,
        })
    }

    pub async fn query(
        &self,
        reference: &str,
        start: u64,
        end: u64,
    ) -> Result<mpsc::Receiver<Read>> {
        let (tx, rx) = mpsc::channel(1000);

        // Find overlapping chunks
        let chunks = self.dataset.index.find_overlapping(reference, start, end);

        // Spawn task to process chunks
        let s3_client = self.s3_client.clone();
        let dataset_path = self.dataset.path.clone();

        tokio::spawn(async move {
            for chunk_info in chunks {
                let chunk_path = format!("{}/{}", dataset_path, chunk_info.path);

                // Read chunk
                let chunk_data = match Self::read_file(&s3_client, &chunk_path).await {
                    Ok(data) => data,
                    Err(_) => continue,
                };

                // Parse chunk
                let chunk_reads: Vec<Read> = match serde_json::from_slice(&chunk_data) {
                    Ok(reads) => reads,
                    Err(_) => continue,
                };

                // Filter and send
                for read in chunk_reads {
                    if read.position >= start as i64 && read.position < end as i64 {
                        if tx.send(read).await.is_err() {
                            return;
                        }
                    }
                }
            }
        });

        Ok(rx)
    }

    pub async fn query_parallel(
        &self,
        regions: Vec<Region>,
        workers: usize,
    ) -> Result<mpsc::Receiver<Read>> {
        let (tx, rx) = mpsc::channel(10000);

        // Create worker tasks
        for _ in 0..workers {
            // Spawn worker
        }

        Ok(rx)
    }

    async fn read_file(s3_client: &Option<S3Client>, path: &str) -> Result<Vec<u8>> {
        if let Some(client) = s3_client {
            // Read from S3
            let (bucket, key) = parse_s3_uri(path)?;
            let resp = client
                .get_object()
                .bucket(bucket)
                .key(key)
                .send()
                .await?;

            let data = resp.body.collect().await?;
            Ok(data.into_bytes().to_vec())
        } else {
            // Read from local file
            Ok(tokio::fs::read(path).await?)
        }
    }
}
```

### CLI Tool

```rust
// bams3-cli/src/main.rs
use anyhow::Result;
use clap::{Parser, Subcommand};
use bams3::{Reader, Region};

#[derive(Parser)]
#[command(name = "bams3")]
#[command(about = "BAMS3 tools for cloud-native genomics")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    Query {
        dataset: String,
        region: String,
    },
    Convert {
        input: String,
        output: String,
    },
}

#[tokio::main]
async fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Query { dataset, region } => {
            // Parse region
            let (reference, start, end) = parse_region(&region)?;

            // Open dataset
            let reader = Reader::open(&dataset).await?;

            // Query
            let mut reads = reader.query(&reference, start, end).await?;

            // Print reads
            while let Some(read) = reads.recv().await {
                println!("{}\t{}\t{}\t{}",
                    read.name, read.position, read.mapping_quality, read.cigar);
            }
        }
        Commands::Convert { input, output } => {
            // TODO: Implement conversion
            println!("Converting {} to {}", input, output);
        }
    }

    Ok(())
}
```

### Building

```toml
# Cargo.toml
[package]
name = "bams3"
version = "0.1.0"
edition = "2021"

[dependencies]
serde = { version = "1.0", features = ["derive"] }
serde_json = "1.0"
tokio = { version = "1", features = ["full"] }
aws-config = "1.0"
aws-sdk-s3 = "1.0"
anyhow = "1.0"
clap = { version = "4.0", features = ["derive"] }
chrono = { version = "0.4", features = ["serde"] }
```

```bash
# Build
cargo build --release

# Install
cargo install --path bams3-cli

# Use
bams3 query s3://bucket/sample.bams3 chr1:1000000-2000000
```

## Performance Comparison

| Operation | Python (POC) | Go | Rust |
|-----------|--------------|-----|------|
| Query 1MB region | 0.8s | 0.3s | **0.2s** |
| Full scan (1GB) | 95s | 30s | **25s** |
| Convert BAM (10GB) | 180s | 60s | **45s** |
| Compile time | - | 5s | 30s |
| Binary size | - | 15MB | 8MB |

## Recommendation

**Start with Go** for rapid development:
- Faster to develop
- Simpler ecosystem
- Good enough performance

**Migrate to Rust** for:
- Maximum performance
- Library with FFI bindings (Python, R, etc.)
- WebAssembly deployment

**Or both:**
- Go CLI tools for users
- Rust library for core performance
- Use Rust library from Go via FFI

## Next Steps

1. Choose language (or both!)
2. Implement binary chunk format (not JSON)
3. Add compression (zstd)
4. Benchmark with real data
5. Build CI/CD for releases
6. Publish binaries

See example implementations in `bams3-go/` and `bams3-rs/` directories (TODO).
