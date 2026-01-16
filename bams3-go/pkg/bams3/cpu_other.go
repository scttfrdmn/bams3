//go:build !darwin && !linux

package bams3

import "runtime"

// detectOptimalWorkers fallback for unsupported operating systems
func detectOptimalWorkers() int {
	return runtime.NumCPU()
}
