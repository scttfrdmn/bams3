#!/usr/bin/env python3
"""
Benchmark different S3 access methods for genomics workflows.

Compares:
1. Copy-then-process (baseline)
2. FUSE mounting (mountpoint-s3, goofys, s3fs)
3. Streaming
4. Direct S3 access (if tool supports)

Usage:
    python benchmark.py --config benchmark-config.yaml
    python benchmark.py --s3-file s3://bucket/test.bam --operation flagstat
"""

import argparse
import subprocess
import time
import tempfile
import shutil
import os
import sys
import json
from pathlib import Path
from dataclasses import dataclass, asdict
from typing import List, Dict, Optional
import yaml


@dataclass
class BenchmarkResult:
    """Results from a single benchmark run."""
    method: str
    operation: str
    file_size_bytes: int
    duration_seconds: float
    success: bool
    error_message: Optional[str] = None
    throughput_mbps: Optional[float] = None
    peak_memory_mb: Optional[float] = None

    def __post_init__(self):
        if self.success and self.duration_seconds > 0:
            self.throughput_mbps = (self.file_size_bytes / (1024 * 1024)) / self.duration_seconds


class Benchmarker:
    """Run benchmarks comparing different S3 access methods."""

    def __init__(self, s3_file: str, local_cache_dir: str = "/tmp/s3-benchmark"):
        self.s3_file = s3_file
        self.local_cache_dir = Path(local_cache_dir)
        self.local_cache_dir.mkdir(parents=True, exist_ok=True)
        self.file_size = self._get_file_size()
        self.results: List[BenchmarkResult] = []

    def _get_file_size(self) -> int:
        """Get S3 file size in bytes."""
        # Parse s3://bucket/key
        s3_path = self.s3_file[5:]  # Remove s3://
        bucket, key = s3_path.split('/', 1)

        try:
            result = subprocess.run(
                ['aws', 's3api', 'head-object', '--bucket', bucket, '--key', key],
                capture_output=True,
                text=True,
                check=True
            )
            import json
            metadata = json.loads(result.stdout)
            return metadata['ContentLength']
        except Exception as e:
            print(f"Warning: Could not get file size: {e}")
            return 0

    def _run_command(self, cmd: List[str], timeout: int = 600) -> tuple[bool, float, Optional[str]]:
        """
        Run a command and return (success, duration, error_message).

        Args:
            cmd: Command to run
            timeout: Timeout in seconds

        Returns:
            (success, duration, error_message)
        """
        start = time.time()
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout,
                check=True
            )
            duration = time.time() - start
            return True, duration, None
        except subprocess.TimeoutExpired:
            duration = time.time() - start
            return False, duration, f"Timeout after {timeout}s"
        except subprocess.CalledProcessError as e:
            duration = time.time() - start
            return False, duration, f"Command failed: {e.stderr}"
        except Exception as e:
            duration = time.time() - start
            return False, duration, str(e)

    def benchmark_copy_then_process(self, operation: str) -> BenchmarkResult:
        """
        Benchmark 1: Copy file locally, then process.
        This is the baseline that many users do today.
        """
        print(f"  Running copy-then-process benchmark...")

        with tempfile.TemporaryDirectory(dir=self.local_cache_dir) as tmpdir:
            local_file = Path(tmpdir) / "data_file"

            # Copy from S3
            copy_cmd = ['aws', 's3', 'cp', self.s3_file, str(local_file)]
            success, copy_time, error = self._run_command(copy_cmd)

            if not success:
                return BenchmarkResult(
                    method="copy-then-process",
                    operation=operation,
                    file_size_bytes=self.file_size,
                    duration_seconds=copy_time,
                    success=False,
                    error_message=f"Copy failed: {error}"
                )

            # Process locally
            process_cmd = self._get_operation_command(operation, str(local_file))
            success, process_time, error = self._run_command(process_cmd)

            total_time = copy_time + process_time

            return BenchmarkResult(
                method="copy-then-process",
                operation=operation,
                file_size_bytes=self.file_size,
                duration_seconds=total_time,
                success=success,
                error_message=error
            )

    def benchmark_fuse_mount(self, operation: str, mount_tool: str = "mountpoint-s3") -> BenchmarkResult:
        """
        Benchmark 2: FUSE mount and process.

        Args:
            operation: Operation to perform
            mount_tool: FUSE tool to use (mountpoint-s3, goofys, s3fs)
        """
        print(f"  Running FUSE mount benchmark ({mount_tool})...")

        # Parse S3 path
        s3_path = self.s3_file[5:]  # Remove s3://
        bucket, key = s3_path.split('/', 1)

        with tempfile.TemporaryDirectory(dir=self.local_cache_dir) as tmpdir:
            mount_point = Path(tmpdir) / "mount"
            mount_point.mkdir()
            cache_dir = Path(tmpdir) / "cache"
            cache_dir.mkdir()

            # Mount command depends on tool
            if mount_tool == "mountpoint-s3":
                mount_cmd = [
                    'mount-s3',
                    bucket,
                    str(mount_point),
                    '--cache', str(cache_dir),
                    '--read-only'
                ]
            elif mount_tool == "goofys":
                mount_cmd = ['goofys', bucket, str(mount_point)]
            elif mount_tool == "s3fs":
                mount_cmd = [
                    's3fs',
                    bucket,
                    str(mount_point),
                    '-o', f'use_cache={cache_dir}',
                    '-o', 'ro'
                ]
            else:
                return BenchmarkResult(
                    method=f"fuse-{mount_tool}",
                    operation=operation,
                    file_size_bytes=self.file_size,
                    duration_seconds=0,
                    success=False,
                    error_message=f"Unknown mount tool: {mount_tool}"
                )

            # Mount
            try:
                subprocess.run(mount_cmd, check=True, capture_output=True, timeout=30)
            except Exception as e:
                return BenchmarkResult(
                    method=f"fuse-{mount_tool}",
                    operation=operation,
                    file_size_bytes=self.file_size,
                    duration_seconds=0,
                    success=False,
                    error_message=f"Mount failed: {e}"
                )

            try:
                # Process through mount
                mounted_file = mount_point / key
                process_cmd = self._get_operation_command(operation, str(mounted_file))
                success, duration, error = self._run_command(process_cmd)

                return BenchmarkResult(
                    method=f"fuse-{mount_tool}",
                    operation=operation,
                    file_size_bytes=self.file_size,
                    duration_seconds=duration,
                    success=success,
                    error_message=error
                )
            finally:
                # Unmount
                try:
                    subprocess.run(['umount', str(mount_point)], timeout=10)
                except:
                    # Try fusermount if umount fails
                    try:
                        subprocess.run(['fusermount', '-u', str(mount_point)], timeout=10)
                    except:
                        pass

    def benchmark_streaming(self, operation: str) -> BenchmarkResult:
        """
        Benchmark 3: Stream from S3 (if operation supports it).
        """
        print(f"  Running streaming benchmark...")

        # Only some operations support streaming
        if not self._operation_supports_streaming(operation):
            return BenchmarkResult(
                method="streaming",
                operation=operation,
                file_size_bytes=self.file_size,
                duration_seconds=0,
                success=False,
                error_message="Operation does not support streaming"
            )

        # Build streaming command
        stream_cmd = self._get_streaming_command(operation)
        success, duration, error = self._run_command(stream_cmd)

        return BenchmarkResult(
            method="streaming",
            operation=operation,
            file_size_bytes=self.file_size,
            duration_seconds=duration,
            success=success,
            error_message=error
        )

    def benchmark_direct_s3(self, operation: str) -> BenchmarkResult:
        """
        Benchmark 4: Direct S3 access (if tool supports s3:// URIs).
        """
        print(f"  Running direct S3 benchmark...")

        # Check if tool supports S3
        if not self._tool_supports_s3(operation):
            return BenchmarkResult(
                method="direct-s3",
                operation=operation,
                file_size_bytes=self.file_size,
                duration_seconds=0,
                success=False,
                error_message="Tool does not support s3:// URIs"
            )

        # Process with s3:// URI
        process_cmd = self._get_operation_command(operation, self.s3_file)
        success, duration, error = self._run_command(process_cmd)

        return BenchmarkResult(
            method="direct-s3",
            operation=operation,
            file_size_bytes=self.file_size,
            duration_seconds=duration,
            success=success,
            error_message=error
        )

    def _get_operation_command(self, operation: str, file_path: str) -> List[str]:
        """Get the command for a specific operation."""
        # Map operation names to actual commands
        commands = {
            'flagstat': ['samtools', 'flagstat', file_path],
            'view-head': ['samtools', 'view', file_path],  # Will pipe to head
            'idxstats': ['samtools', 'idxstats', file_path],
            'count-reads': ['samtools', 'view', '-c', file_path],
            'region-query': ['samtools', 'view', file_path, 'chr1:1000000-2000000'],
            'head': ['head', '-n', '1000', file_path],
            'wc': ['wc', '-l', file_path],
        }

        cmd = commands.get(operation)
        if cmd is None:
            raise ValueError(f"Unknown operation: {operation}")

        # Redirect output to /dev/null to avoid I/O overhead
        return cmd + ['>', '/dev/null'] if '>' not in cmd else cmd

    def _get_streaming_command(self, operation: str) -> List[str]:
        """Get streaming command for operation."""
        base = ['aws', 's3', 'cp', self.s3_file, '-']

        if operation == 'flagstat':
            return base + ['|', 'samtools', 'flagstat', '-']
        elif operation == 'count-reads':
            return base + ['|', 'samtools', 'view', '-c', '-']
        elif operation == 'head':
            return base + ['|', 'head', '-n', '1000']
        elif operation == 'wc':
            return base + ['|', 'wc', '-l']
        else:
            raise ValueError(f"Streaming not supported for: {operation}")

    def _operation_supports_streaming(self, operation: str) -> bool:
        """Check if operation can be done via streaming."""
        streaming_ops = {'flagstat', 'count-reads', 'head', 'wc'}
        return operation in streaming_ops

    def _tool_supports_s3(self, operation: str) -> bool:
        """Check if the tool supports S3 URIs."""
        # samtools/bcftools with htslib compiled with S3 support
        samtools_ops = {'flagstat', 'view-head', 'idxstats', 'count-reads', 'region-query'}

        if operation in samtools_ops:
            # Check if samtools has S3 support
            try:
                result = subprocess.run(
                    ['samtools', 'version'],
                    capture_output=True,
                    text=True
                )
                return 's3' in result.stdout.lower()
            except:
                return False

        return False

    def run_all_benchmarks(self, operation: str) -> List[BenchmarkResult]:
        """Run all applicable benchmarks for an operation."""
        print(f"\nBenchmarking operation: {operation}")
        print(f"File: {self.s3_file}")
        print(f"Size: {self.file_size / (1024**3):.2f} GB\n")

        results = []

        # Always run copy-then-process (baseline)
        results.append(self.benchmark_copy_then_process(operation))

        # Try FUSE mounts (check which are available)
        for mount_tool in ['mountpoint-s3', 'goofys', 's3fs']:
            if shutil.which(mount_tool) or (mount_tool == 'mountpoint-s3' and shutil.which('mount-s3')):
                results.append(self.benchmark_fuse_mount(operation, mount_tool))
            else:
                print(f"  Skipping {mount_tool} (not installed)")

        # Try streaming
        if self._operation_supports_streaming(operation):
            results.append(self.benchmark_streaming(operation))
        else:
            print(f"  Skipping streaming (not supported for {operation})")

        # Try direct S3
        if self._tool_supports_s3(operation):
            results.append(self.benchmark_direct_s3(operation))
        else:
            print(f"  Skipping direct S3 (tool doesn't support s3:// URIs)")

        self.results.extend(results)
        return results

    def print_results(self):
        """Print results in a readable table format."""
        if not self.results:
            print("No results to display")
            return

        print("\n" + "="*80)
        print("BENCHMARK RESULTS")
        print("="*80)
        print(f"\nFile: {self.s3_file}")
        print(f"Size: {self.file_size / (1024**3):.2f} GB")
        print()

        # Group by operation
        operations = {}
        for result in self.results:
            if result.operation not in operations:
                operations[result.operation] = []
            operations[result.operation].append(result)

        for operation, results in operations.items():
            print(f"\nOperation: {operation}")
            print("-" * 80)
            print(f"{'Method':<25} {'Duration':<12} {'Throughput':<15} {'Status':<10}")
            print("-" * 80)

            # Sort by duration (successful ones first)
            results_sorted = sorted(results, key=lambda r: (not r.success, r.duration_seconds))

            for result in results_sorted:
                duration_str = f"{result.duration_seconds:.2f}s" if result.success else "FAILED"
                throughput_str = f"{result.throughput_mbps:.1f} MB/s" if result.throughput_mbps else "N/A"
                status = "✓" if result.success else f"✗ {result.error_message}"

                print(f"{result.method:<25} {duration_str:<12} {throughput_str:<15} {status:<10}")

            # Calculate speedup vs baseline
            baseline = next((r for r in results if r.method == "copy-then-process" and r.success), None)
            if baseline:
                print()
                print("Speedup vs copy-then-process:")
                for result in results_sorted:
                    if result.success and result.method != "copy-then-process":
                        speedup = baseline.duration_seconds / result.duration_seconds
                        print(f"  {result.method:<25} {speedup:.2f}x")

        print("\n" + "="*80)

    def save_results(self, output_file: str):
        """Save results to JSON file."""
        data = {
            's3_file': self.s3_file,
            'file_size_bytes': self.file_size,
            'results': [asdict(r) for r in self.results]
        }

        with open(output_file, 'w') as f:
            json.dump(data, f, indent=2)

        print(f"\nResults saved to: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Benchmark S3 access methods for genomics workflows'
    )
    parser.add_argument(
        '--s3-file',
        required=True,
        help='S3 file to benchmark (e.g., s3://bucket/file.bam)'
    )
    parser.add_argument(
        '--operation',
        default='flagstat',
        choices=['flagstat', 'view-head', 'idxstats', 'count-reads', 'region-query', 'head', 'wc'],
        help='Operation to benchmark'
    )
    parser.add_argument(
        '--cache-dir',
        default='/tmp/s3-benchmark',
        help='Directory for temporary files and cache'
    )
    parser.add_argument(
        '--output',
        default='benchmark-results.json',
        help='Output file for results (JSON)'
    )

    args = parser.parse_args()

    # Verify S3 file exists
    if not args.s3_file.startswith('s3://'):
        print("Error: --s3-file must be an S3 URI (s3://bucket/key)")
        sys.exit(1)

    benchmarker = Benchmarker(args.s3_file, args.cache_dir)
    benchmarker.run_all_benchmarks(args.operation)
    benchmarker.print_results()
    benchmarker.save_results(args.output)


if __name__ == '__main__':
    main()
