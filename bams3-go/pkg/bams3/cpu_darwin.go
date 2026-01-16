//go:build darwin

package bams3

import (
	"runtime"
	"syscall"
)

// detectOptimalWorkers detects optimal worker count on macOS
// Prefers performance cores over efficiency cores for CPU-intensive tasks
func detectOptimalWorkers() int {
	// Try to detect performance cores (Apple Silicon M1/M2/M3)
	perfCores := getPerfCores()
	if perfCores > 0 {
		return perfCores
	}

	// Fallback to all logical CPUs
	return runtime.NumCPU()
}

// getPerfCores returns the number of performance cores on Apple Silicon
func getPerfCores() int {
	// Try hw.perflevel0.physicalcpu (performance cores)
	// Note: syscall.Sysctl returns raw bytes, not a string
	result, err := syscall.Sysctl("hw.perflevel0.physicalcpu")
	if err == nil && len(result) > 0 {
		// Parse as little-endian integer
		count := int(result[0])
		if len(result) > 1 {
			count |= int(result[1]) << 8
		}
		if count > 0 {
			return count
		}
	}

	// Fallback: try hw.physicalcpu (total physical cores)
	result, err = syscall.Sysctl("hw.physicalcpu")
	if err == nil && len(result) > 0 {
		count := int(result[0])
		if len(result) > 1 {
			count |= int(result[1]) << 8
		}
		if count > 0 {
			return count
		}
	}

	return 0
}
