package bams3

import (
	"context"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/aws/aws-sdk-go-v2/aws"
	"github.com/aws/aws-sdk-go-v2/config"
	"github.com/aws/aws-sdk-go-v2/feature/s3/manager"
	"github.com/aws/aws-sdk-go-v2/service/s3"
)

// Storage is an interface for reading/writing BAMS3 data
// Supports both local filesystem and S3
type Storage interface {
	// Read reads a file
	ReadFile(path string) ([]byte, error)

	// Write writes a file
	WriteFile(path string, data []byte) error

	// List lists files matching a prefix
	List(prefix string) ([]string, error)

	// Exists checks if a file exists
	Exists(path string) (bool, error)

	// MkdirAll creates directory structure
	MkdirAll(path string) error

	// GetBasePath returns the base path
	GetBasePath() string

	// IsS3 returns true if this is S3 storage
	IsS3() bool
}

// LocalStorage implements Storage for local filesystem
type LocalStorage struct {
	basePath string
}

// NewLocalStorage creates a new local storage backend
func NewLocalStorage(basePath string) *LocalStorage {
	return &LocalStorage{basePath: basePath}
}

func (s *LocalStorage) ReadFile(path string) ([]byte, error) {
	fullPath := filepath.Join(s.basePath, path)
	return os.ReadFile(fullPath)
}

func (s *LocalStorage) WriteFile(path string, data []byte) error {
	fullPath := filepath.Join(s.basePath, path)
	// Ensure directory exists
	dir := filepath.Dir(fullPath)
	if err := os.MkdirAll(dir, 0755); err != nil {
		return err
	}
	return os.WriteFile(fullPath, data, 0644)
}

func (s *LocalStorage) List(prefix string) ([]string, error) {
	fullPath := filepath.Join(s.basePath, prefix)
	var files []string

	err := filepath.Walk(fullPath, func(path string, info os.FileInfo, err error) error {
		if err != nil {
			return err
		}
		if !info.IsDir() {
			relPath, err := filepath.Rel(s.basePath, path)
			if err != nil {
				return err
			}
			files = append(files, relPath)
		}
		return nil
	})

	return files, err
}

func (s *LocalStorage) Exists(path string) (bool, error) {
	fullPath := filepath.Join(s.basePath, path)
	_, err := os.Stat(fullPath)
	if err == nil {
		return true, nil
	}
	if os.IsNotExist(err) {
		return false, nil
	}
	return false, err
}

func (s *LocalStorage) MkdirAll(path string) error {
	fullPath := filepath.Join(s.basePath, path)
	return os.MkdirAll(fullPath, 0755)
}

func (s *LocalStorage) GetBasePath() string {
	return s.basePath
}

func (s *LocalStorage) IsS3() bool {
	return false
}

// S3Storage implements Storage for AWS S3
type S3Storage struct {
	bucket     string
	prefix     string
	client     *s3.Client
	uploader   *manager.Uploader
	downloader *manager.Downloader
	ctx        context.Context
}

// NewS3Storage creates a new S3 storage backend
// path should be in format: s3://bucket/prefix
func NewS3Storage(path string) (*S3Storage, error) {
	// Parse S3 path
	if !strings.HasPrefix(path, "s3://") {
		return nil, fmt.Errorf("invalid S3 path: %s (must start with s3://)", path)
	}

	path = strings.TrimPrefix(path, "s3://")
	parts := strings.SplitN(path, "/", 2)
	bucket := parts[0]
	prefix := ""
	if len(parts) > 1 {
		prefix = parts[1]
	}

	// Load AWS config
	ctx := context.Background()
	cfg, err := config.LoadDefaultConfig(ctx)
	if err != nil {
		return nil, fmt.Errorf("failed to load AWS config: %w", err)
	}

	client := s3.NewFromConfig(cfg)

	return &S3Storage{
		bucket:     bucket,
		prefix:     prefix,
		client:     client,
		uploader:   manager.NewUploader(client),
		downloader: manager.NewDownloader(client),
		ctx:        ctx,
	}, nil
}

func (s *S3Storage) getFullKey(path string) string {
	if s.prefix == "" {
		return path
	}
	return s.prefix + "/" + path
}

func (s *S3Storage) ReadFile(path string) ([]byte, error) {
	key := s.getFullKey(path)

	// Download to memory
	buf := manager.NewWriteAtBuffer([]byte{})
	_, err := s.downloader.Download(s.ctx, buf, &s3.GetObjectInput{
		Bucket: aws.String(s.bucket),
		Key:    aws.String(key),
	})
	if err != nil {
		return nil, fmt.Errorf("failed to download s3://%s/%s: %w", s.bucket, key, err)
	}

	return buf.Bytes(), nil
}

func (s *S3Storage) WriteFile(path string, data []byte) error {
	key := s.getFullKey(path)

	_, err := s.uploader.Upload(s.ctx, &s3.PutObjectInput{
		Bucket: aws.String(s.bucket),
		Key:    aws.String(key),
		Body:   strings.NewReader(string(data)),
	})
	if err != nil {
		return fmt.Errorf("failed to upload to s3://%s/%s: %w", s.bucket, key, err)
	}

	return nil
}

func (s *S3Storage) List(prefix string) ([]string, error) {
	fullPrefix := s.getFullKey(prefix)

	var files []string
	paginator := s3.NewListObjectsV2Paginator(s.client, &s3.ListObjectsV2Input{
		Bucket: aws.String(s.bucket),
		Prefix: aws.String(fullPrefix),
	})

	for paginator.HasMorePages() {
		page, err := paginator.NextPage(s.ctx)
		if err != nil {
			return nil, fmt.Errorf("failed to list objects: %w", err)
		}

		for _, obj := range page.Contents {
			key := aws.ToString(obj.Key)
			// Remove the full prefix to get relative path
			if s.prefix != "" {
				key = strings.TrimPrefix(key, s.prefix+"/")
			}
			files = append(files, key)
		}
	}

	return files, nil
}

func (s *S3Storage) Exists(path string) (bool, error) {
	key := s.getFullKey(path)

	_, err := s.client.HeadObject(s.ctx, &s3.HeadObjectInput{
		Bucket: aws.String(s.bucket),
		Key:    aws.String(key),
	})
	if err != nil {
		// Check if it's a not found error
		if strings.Contains(err.Error(), "NotFound") || strings.Contains(err.Error(), "404") {
			return false, nil
		}
		return false, err
	}

	return true, nil
}

func (s *S3Storage) MkdirAll(path string) error {
	// S3 doesn't have directories, so this is a no-op
	return nil
}

func (s *S3Storage) GetBasePath() string {
	if s.prefix == "" {
		return fmt.Sprintf("s3://%s", s.bucket)
	}
	return fmt.Sprintf("s3://%s/%s", s.bucket, s.prefix)
}

func (s *S3Storage) IsS3() bool {
	return true
}

// NewStorage creates the appropriate storage backend based on path
func NewStorage(path string) (Storage, error) {
	if strings.HasPrefix(path, "s3://") {
		return NewS3Storage(path)
	}
	return NewLocalStorage(path), nil
}

// FileInfo represents file metadata
type FileInfo struct {
	Path    string
	Size    int64
	ModTime time.Time
}

// StreamReader provides streaming access to files
type StreamReader interface {
	io.ReadCloser
}

// StreamWriter provides streaming write access
type StreamWriter interface {
	io.WriteCloser
}
