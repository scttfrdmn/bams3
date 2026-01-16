package bams3

import (
	"bytes"
	"context"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"
	"sync"

	"github.com/aws/aws-sdk-go-v2/aws"
	"github.com/aws/aws-sdk-go-v2/config"
	"github.com/aws/aws-sdk-go-v2/feature/s3/manager"
	"github.com/aws/aws-sdk-go-v2/service/s3"
)

// S3Writer handles uploading BAMS3 format to S3
type S3Writer struct {
	client   *s3.Client
	uploader *manager.Uploader
	bucket   string
	prefix   string
	region   string

	// Progress tracking
	uploadedBytes int64
	uploadMutex   sync.Mutex
}

// S3URI represents a parsed S3 URI
type S3URI struct {
	Bucket string
	Prefix string
}

// ParseS3URI parses an S3 URI like s3://bucket/path/to/object
func ParseS3URI(uri string) (*S3URI, error) {
	if !strings.HasPrefix(uri, "s3://") {
		return nil, fmt.Errorf("invalid S3 URI: must start with s3://")
	}

	// Remove s3:// prefix
	path := strings.TrimPrefix(uri, "s3://")

	// Split into bucket and prefix
	parts := strings.SplitN(path, "/", 2)
	if len(parts) == 0 || parts[0] == "" {
		return nil, fmt.Errorf("invalid S3 URI: missing bucket name")
	}

	bucket := parts[0]
	prefix := ""
	if len(parts) == 2 {
		prefix = parts[1]
	}

	return &S3URI{
		Bucket: bucket,
		Prefix: prefix,
	}, nil
}

// IsS3URI checks if a path is an S3 URI
func IsS3URI(path string) bool {
	return strings.HasPrefix(path, "s3://")
}

// NewS3Writer creates a new S3 writer
func NewS3Writer(ctx context.Context, s3URI string, region string) (*S3Writer, error) {
	// Parse S3 URI
	uri, err := ParseS3URI(s3URI)
	if err != nil {
		return nil, err
	}

	// Load AWS configuration
	cfg, err := config.LoadDefaultConfig(ctx,
		config.WithRegion(region),
	)
	if err != nil {
		return nil, fmt.Errorf("failed to load AWS config: %w", err)
	}

	// If region not specified, try to detect from bucket
	if region == "" {
		region = cfg.Region
	}

	// Create S3 client
	client := s3.NewFromConfig(cfg)

	// Create uploader with custom settings
	uploader := manager.NewUploader(client, func(u *manager.Uploader) {
		// Use 10MB part size (optimal for most workloads)
		u.PartSize = 10 * 1024 * 1024
		// Allow up to 3 concurrent uploads
		u.Concurrency = 3
	})

	return &S3Writer{
		client:   client,
		uploader: uploader,
		bucket:   uri.Bucket,
		prefix:   uri.Prefix,
		region:   region,
	}, nil
}

// WriteFile uploads a file to S3
func (w *S3Writer) WriteFile(ctx context.Context, relativePath string, data []byte) error {
	// Construct full S3 key
	key := filepath.Join(w.prefix, relativePath)
	// Normalize path separators for S3 (always use /)
	key = strings.ReplaceAll(key, "\\", "/")

	// Upload using multipart uploader
	reader := bytes.NewReader(data)

	input := &s3.PutObjectInput{
		Bucket: aws.String(w.bucket),
		Key:    aws.String(key),
		Body:   reader,
	}

	_, err := w.uploader.Upload(ctx, input)
	if err != nil {
		return fmt.Errorf("failed to upload %s: %w", key, err)
	}

	// Track uploaded bytes
	w.uploadMutex.Lock()
	w.uploadedBytes += int64(len(data))
	w.uploadMutex.Unlock()

	return nil
}

// WriteFileStream uploads a file from a reader
func (w *S3Writer) WriteFileStream(ctx context.Context, relativePath string, reader io.Reader, size int64) error {
	// Construct full S3 key
	key := filepath.Join(w.prefix, relativePath)
	key = strings.ReplaceAll(key, "\\", "/")

	input := &s3.PutObjectInput{
		Bucket: aws.String(w.bucket),
		Key:    aws.String(key),
		Body:   reader,
	}

	_, err := w.uploader.Upload(ctx, input)
	if err != nil {
		return fmt.Errorf("failed to upload %s: %w", key, err)
	}

	// Track uploaded bytes
	w.uploadMutex.Lock()
	w.uploadedBytes += size
	w.uploadMutex.Unlock()

	return nil
}

// GetUploadedBytes returns total bytes uploaded
func (w *S3Writer) GetUploadedBytes() int64 {
	w.uploadMutex.Lock()
	defer w.uploadMutex.Unlock()
	return w.uploadedBytes
}

// S3ChunkWriter adapts S3Writer to work with ParallelWriter
type S3ChunkWriter struct {
	s3Writer *S3Writer
	ctx      context.Context
}

// NewS3ChunkWriter creates a chunk writer for S3
func NewS3ChunkWriter(ctx context.Context, s3URI string, region string) (*S3ChunkWriter, error) {
	writer, err := NewS3Writer(ctx, s3URI, region)
	if err != nil {
		return nil, err
	}

	return &S3ChunkWriter{
		s3Writer: writer,
		ctx:      ctx,
	}, nil
}

// WriteChunk writes a chunk to S3
func (w *S3ChunkWriter) WriteChunk(chunkPath string, data []byte) error {
	return w.s3Writer.WriteFile(w.ctx, chunkPath, data)
}

// WriteMetadata writes metadata file to S3
func (w *S3ChunkWriter) WriteMetadata(filename string, data []byte) error {
	return w.s3Writer.WriteFile(w.ctx, filename, data)
}

// Close finalizes the S3 writer
func (w *S3ChunkWriter) Close() error {
	// Log final upload stats
	uploaded := w.s3Writer.GetUploadedBytes()
	fmt.Fprintf(os.Stderr, "\nTotal uploaded to S3: %.1f MB\n", float64(uploaded)/(1024*1024))
	return nil
}

// GetS3URI returns the S3 URI for this writer
func (w *S3ChunkWriter) GetS3URI() string {
	if w.s3Writer.prefix != "" {
		return fmt.Sprintf("s3://%s/%s", w.s3Writer.bucket, w.s3Writer.prefix)
	}
	return fmt.Sprintf("s3://%s", w.s3Writer.bucket)
}
