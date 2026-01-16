//go:build linux

package bams3

import (
	"bufio"
	"os"
	"runtime"
	"strconv"
	"strings"
)

// detectOptimalWorkers detects optimal worker count on Linux
// Attempts to detect performance cores vs efficiency cores
func detectOptimalWorkers() int {
	// Try to detect Intel hybrid architecture (P-cores vs E-cores)
	// by looking at CPU frequencies
	perfCores := detectPerfCoresLinux()
	if perfCores > 0 {
		return perfCores
	}

	// Fallback to all logical CPUs
	return runtime.NumCPU()
}

// detectPerfCoresLinux attempts to detect performance cores on Linux
func detectPerfCoresLinux() int {
	// Read /proc/cpuinfo to analyze core types
	file, err := os.Open("/proc/cpuinfo")
	if err != nil {
		return 0
	}
	defer file.Close()

	type cpuCore struct {
		id        int
		coreID    int
		frequency float64
	}

	cores := make(map[int]cpuCore)
	currentCPU := -1
	currentCoreID := -1

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		fields := strings.SplitN(line, ":", 2)
		if len(fields) != 2 {
			continue
		}

		key := strings.TrimSpace(fields[0])
		value := strings.TrimSpace(fields[1])

		switch key {
		case "processor":
			if id, err := strconv.Atoi(value); err == nil {
				currentCPU = id
			}
		case "core id":
			if id, err := strconv.Atoi(value); err == nil {
				currentCoreID = id
			}
		case "cpu MHz":
			if freq, err := strconv.ParseFloat(value, 64); err == nil && currentCPU >= 0 {
				cores[currentCPU] = cpuCore{
					id:        currentCPU,
					coreID:    currentCoreID,
					frequency: freq,
				}
			}
		}
	}

	// If we have frequency data, identify high-performance cores
	// P-cores typically run at higher base frequencies
	if len(cores) > 0 {
		// Find unique core IDs with their max frequency
		uniqueCores := make(map[int]float64)
		for _, core := range cores {
			if core.coreID >= 0 {
				if freq, exists := uniqueCores[core.coreID]; !exists || core.frequency > freq {
					uniqueCores[core.coreID] = core.frequency
				}
			}
		}

		// If we have variation in frequencies, use the high-frequency cores
		// For homogeneous systems, this will just return all cores
		if len(uniqueCores) > 0 {
			// Calculate median frequency
			var freqs []float64
			for _, freq := range uniqueCores {
				freqs = append(freqs, freq)
			}

			// Count cores above median (likely P-cores)
			if len(freqs) > 2 {
				var sum float64
				for _, f := range freqs {
					sum += f
				}
				avgFreq := sum / float64(len(freqs))

				perfCount := 0
				for _, freq := range uniqueCores {
					if freq >= avgFreq*0.9 { // Within 10% of average
						perfCount++
					}
				}

				if perfCount > 0 && perfCount < len(uniqueCores) {
					// We detected a hybrid architecture
					return perfCount
				}
			}
		}
	}

	// No hybrid architecture detected or failed to parse
	return 0
}
