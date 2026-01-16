package bams3

import (
	"bufio"
	"container/heap"
	"encoding/gob"
	"fmt"
	"io"
	"os"
	"path/filepath"
)

// SpillWriter writes sorted reads to a temporary file
type SpillWriter struct {
	file    *os.File
	encoder *gob.Encoder
	writer  *bufio.Writer
	count   int
	bytes   int64
}

// NewSpillWriter creates a new spill writer
func NewSpillWriter(dir string, spillNum int) (*SpillWriter, error) {
	path := filepath.Join(dir, fmt.Sprintf("spill-%04d.dat", spillNum))
	file, err := os.Create(path)
	if err != nil {
		return nil, fmt.Errorf("failed to create spill file: %w", err)
	}

	writer := bufio.NewWriterSize(file, 1*MB)
	encoder := gob.NewEncoder(writer)

	return &SpillWriter{
		file:    file,
		encoder: encoder,
		writer:  writer,
		count:   0,
		bytes:   0,
	}, nil
}

// WriteRead writes a single read to the spill file
func (sw *SpillWriter) WriteRead(read *Read) error {
	if err := sw.encoder.Encode(read); err != nil {
		return fmt.Errorf("failed to encode read: %w", err)
	}
	sw.count++
	sw.bytes += int64(len(read.Sequence) + len(read.Quality) + 100)
	return nil
}

// Close finalizes and closes the spill file
func (sw *SpillWriter) Close() (SpillFile, error) {
	if err := sw.writer.Flush(); err != nil {
		sw.file.Close()
		return SpillFile{}, fmt.Errorf("failed to flush writer: %w", err)
	}

	path := sw.file.Name()
	if err := sw.file.Close(); err != nil {
		return SpillFile{}, fmt.Errorf("failed to close file: %w", err)
	}

	// Get file size
	stat, err := os.Stat(path)
	if err != nil {
		return SpillFile{}, fmt.Errorf("failed to stat file: %w", err)
	}

	return SpillFile{
		Path:      path,
		NumReads:  sw.count,
		SizeBytes: stat.Size(),
	}, nil
}

// SpillReader reads sorted reads from a spill file
type SpillReader struct {
	file    *os.File
	decoder *gob.Decoder
	reader  *bufio.Reader
	current *Read
	err     error
	eof     bool
}

// NewSpillReader creates a new spill reader
func NewSpillReader(path string) (*SpillReader, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("failed to open spill file: %w", err)
	}

	reader := bufio.NewReaderSize(file, 1*MB)
	decoder := gob.NewDecoder(reader)

	sr := &SpillReader{
		file:    file,
		decoder: decoder,
		reader:  reader,
	}

	// Prime the reader by loading first record
	sr.advance()

	return sr, nil
}

// advance reads the next read from the file
func (sr *SpillReader) advance() {
	if sr.eof {
		return
	}

	var read Read
	err := sr.decoder.Decode(&read)
	if err != nil {
		if err == io.EOF {
			sr.eof = true
			sr.current = nil
		} else {
			sr.err = err
		}
		return
	}

	sr.current = &read
}

// Peek returns the current read without advancing
func (sr *SpillReader) Peek() (*Read, error) {
	if sr.err != nil {
		return nil, sr.err
	}
	if sr.eof {
		return nil, io.EOF
	}
	return sr.current, nil
}

// Next advances to the next read
func (sr *SpillReader) Next() (*Read, error) {
	current := sr.current
	if sr.err != nil {
		return nil, sr.err
	}
	if sr.eof {
		return nil, io.EOF
	}
	sr.advance()
	return current, nil
}

// Close closes the spill reader
func (sr *SpillReader) Close() error {
	return sr.file.Close()
}

// MergeItem represents an item in the k-way merge heap
type MergeItem struct {
	read       *Read
	readerIdx  int
	readerType string // "spill" or "memory"
}

// MergeHeap implements heap.Interface for k-way merge
type MergeHeap []MergeItem

func (h MergeHeap) Len() int { return len(h) }

func (h MergeHeap) Less(i, j int) bool {
	// Sort by reference ID, then position
	if h[i].read.ReferenceID != h[j].read.ReferenceID {
		return h[i].read.ReferenceID < h[j].read.ReferenceID
	}
	return h[i].read.Position < h[j].read.Position
}

func (h MergeHeap) Swap(i, j int) {
	h[i], h[j] = h[j], h[i]
}

func (h *MergeHeap) Push(x interface{}) {
	*h = append(*h, x.(MergeItem))
}

func (h *MergeHeap) Pop() interface{} {
	old := *h
	n := len(old)
	item := old[n-1]
	*h = old[0 : n-1]
	return item
}

// kWayMerge performs a k-way merge of spill files and in-memory buffer
// Returns merged reads in sorted order
func kWayMerge(spillReaders []*SpillReader, memoryBuffer []Read) ([]Read, error) {
	// Initialize heap
	h := &MergeHeap{}
	heap.Init(h)

	// Add first read from each spill file
	for i, reader := range spillReaders {
		read, err := reader.Peek()
		if err != nil && err != io.EOF {
			return nil, fmt.Errorf("failed to peek from spill %d: %w", i, err)
		}
		if err == nil && read != nil {
			heap.Push(h, MergeItem{
				read:       read,
				readerIdx:  i,
				readerType: "spill",
			})
		}
	}

	// Add first read from memory buffer (treated as a separate stream)
	memIdx := 0
	if memIdx < len(memoryBuffer) {
		heap.Push(h, MergeItem{
			read:       &memoryBuffer[memIdx],
			readerIdx:  memIdx,
			readerType: "memory",
		})
	}

	// Merge results
	result := make([]Read, 0, len(memoryBuffer)+1000)

	for h.Len() > 0 {
		// Pop minimum element
		item := heap.Pop(h).(MergeItem)
		result = append(result, *item.read)

		// Refill from the same source
		if item.readerType == "spill" {
			// Advance spill reader
			_, err := spillReaders[item.readerIdx].Next()
			if err != nil && err != io.EOF {
				return nil, fmt.Errorf("failed to read from spill %d: %w", item.readerIdx, err)
			}

			// Add next read if available
			next, err := spillReaders[item.readerIdx].Peek()
			if err != nil && err != io.EOF {
				return nil, fmt.Errorf("failed to peek from spill %d: %w", item.readerIdx, err)
			}
			if err == nil && next != nil {
				heap.Push(h, MergeItem{
					read:       next,
					readerIdx:  item.readerIdx,
					readerType: "spill",
				})
			}
		} else {
			// Advance memory buffer
			memIdx++
			if memIdx < len(memoryBuffer) {
				heap.Push(h, MergeItem{
					read:       &memoryBuffer[memIdx],
					readerIdx:  memIdx,
					readerType: "memory",
				})
			}
		}
	}

	return result, nil
}
