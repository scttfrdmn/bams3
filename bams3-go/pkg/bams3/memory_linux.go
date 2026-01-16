//go:build linux

package bams3

import (
	"bufio"
	"os"
	"strconv"
	"strings"
)

// detectSystemMemory detects actual system memory on Linux
func detectSystemMemory() (total int64, available int64) {
	// Read /proc/meminfo
	file, err := os.Open("/proc/meminfo")
	if err != nil {
		return 0, 0
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	var memTotal, memAvailable int64

	for scanner.Scan() {
		line := scanner.Text()
		fields := strings.Fields(line)
		if len(fields) < 2 {
			continue
		}

		key := strings.TrimSuffix(fields[0], ":")
		value, err := strconv.ParseInt(fields[1], 10, 64)
		if err != nil {
			continue
		}

		// /proc/meminfo reports in KB, convert to bytes
		switch key {
		case "MemTotal":
			memTotal = value * 1024
		case "MemAvailable":
			memAvailable = value * 1024
		}

		// Stop early if we have both values
		if memTotal > 0 && memAvailable > 0 {
			break
		}
	}

	// If MemAvailable is not present (older kernels), estimate from MemFree + Buffers + Cached
	if memTotal > 0 && memAvailable == 0 {
		file.Seek(0, 0)
		scanner = bufio.NewScanner(file)
		var memFree, buffers, cached int64

		for scanner.Scan() {
			line := scanner.Text()
			fields := strings.Fields(line)
			if len(fields) < 2 {
				continue
			}

			key := strings.TrimSuffix(fields[0], ":")
			value, err := strconv.ParseInt(fields[1], 10, 64)
			if err != nil {
				continue
			}

			switch key {
			case "MemFree":
				memFree = value * 1024
			case "Buffers":
				buffers = value * 1024
			case "Cached":
				cached = value * 1024
			}
		}

		memAvailable = memFree + buffers + cached
	}

	return memTotal, memAvailable
}
