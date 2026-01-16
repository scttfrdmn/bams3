//go:build darwin

package bams3

import (
	"syscall"
)

// detectSystemMemory detects actual system memory on macOS
func detectSystemMemory() (total int64, available int64) {
	// Get total physical memory using sysctl
	totalBytes, err := syscall.Sysctl("hw.memsize")
	if err != nil {
		return 0, 0
	}

	// Parse the result as uint64
	var totalMem uint64
	for i := 0; i < len(totalBytes) && i < 8; i++ {
		totalMem |= uint64(totalBytes[i]) << (uint(i) * 8)
	}

	total = int64(totalMem)

	// For available memory, we can use vm_stat or estimate
	// For now, estimate as 75% of total (conservative)
	// In production, this could parse vm_stat output
	available = total * 3 / 4

	return total, available
}
