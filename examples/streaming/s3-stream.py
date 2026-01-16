#!/usr/bin/env python3
"""
Simple utility to stream S3 objects to stdout.

Usage:
    s3-stream.py s3://bucket/key | tool -

This allows POSIX-only tools that accept stdin to process S3 data
without copying the entire file locally first.
"""

import sys
import argparse
import boto3
from botocore.exceptions import ClientError


def parse_s3_uri(uri):
    """Parse s3://bucket/key into components."""
    if not uri.startswith('s3://'):
        raise ValueError(f"Not an S3 URI: {uri}")

    # Remove s3:// prefix
    path = uri[5:]

    # Split into bucket and key
    parts = path.split('/', 1)
    if len(parts) == 1:
        return parts[0], ''
    return parts[0], parts[1]


def stream_s3_to_stdout(s3_uri, chunk_size=1024*1024):
    """
    Stream S3 object to stdout.

    Args:
        s3_uri: S3 URI like s3://bucket/key
        chunk_size: Size of chunks to read (default 1MB)
    """
    bucket, key = parse_s3_uri(s3_uri)

    if not key:
        raise ValueError(f"No key specified in URI: {s3_uri}")

    # Create S3 client
    s3 = boto3.client('s3')

    try:
        # Get object
        response = s3.get_object(Bucket=bucket, Key=key)

        # Stream to stdout
        for chunk in response['Body'].iter_chunks(chunk_size=chunk_size):
            sys.stdout.buffer.write(chunk)

    except ClientError as e:
        error_code = e.response['Error']['Code']
        if error_code == 'NoSuchKey':
            print(f"Error: Object not found: {s3_uri}", file=sys.stderr)
            sys.exit(1)
        elif error_code == 'NoSuchBucket':
            print(f"Error: Bucket not found: {bucket}", file=sys.stderr)
            sys.exit(1)
        else:
            print(f"Error accessing S3: {e}", file=sys.stderr)
            sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description='Stream S3 object to stdout',
        epilog='Example: s3-stream.py s3://bucket/file.fastq.gz | gunzip | tool -'
    )
    parser.add_argument('s3_uri', help='S3 URI (s3://bucket/key)')
    parser.add_argument(
        '--chunk-size',
        type=int,
        default=1024*1024,
        help='Chunk size in bytes (default: 1MB)'
    )

    args = parser.parse_args()

    stream_s3_to_stdout(args.s3_uri, args.chunk_size)


if __name__ == '__main__':
    main()
