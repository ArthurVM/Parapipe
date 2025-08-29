#!/usr/bin/env python3

import sys

def filter_depth(max_depth=5):
    """
    Filters samtools depth output from stdin and outputs to stdout.

    Args:
        min_depth: Minimum read depth to keep.
    """

    for line in sys.stdin:
        chrom, pos, depth = line.strip().split("\t")
        depth = int(depth)
        if depth <= max_depth:
            print(f"{chrom}\t{pos}\t{int(pos)+1}\tN")

if __name__ == "__main__":
    ## Get maximum depth from command line arguments (default: 5)
    if len(sys.argv) > 1:
        try:
            max_depth = int(sys.argv[1])
        except ValueError:
            print("Error: Invalid maximum depth value.", file=sys.stderr)
            sys.exit(1)
    else:
        max_depth = 5

    filter_depth(max_depth)
