//go:build !darwin && !linux

package bams3

// detectSystemMemory fallback for unsupported operating systems
func detectSystemMemory() (total int64, available int64) {
	// Return 0 to trigger fallback to defaults
	return 0, 0
}
