#!/usr/bin/env python3
"""
Exercise 5 - Hamming distance CLI tool.

Requirements implemented:
- Function calculate_hamming_distance(seq1, seq2)
- argparse input: --seq1 and --seq2
- Prints only the final integer to stdout on success
- Handles empty or unequal-length sequences with graceful error reporting
"""

from __future__ import annotations

import argparse
import sys


def calculate_hamming_distance(seq1: str, seq2: str) -> int:
    """Return Hamming distance between two sequences of equal length."""
    if len(seq1) != len(seq2):
        raise ValueError(
            f"Sequences must have the same length (got {len(seq1)} and {len(seq2)})."
        )
    return sum(base1 != base2 for base1, base2 in zip(seq1, seq2))


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Calculate Hamming distance for two DNA sequences."
    )
    parser.add_argument("--seq1", required=True, help="First DNA sequence")
    parser.add_argument("--seq2", required=True, help="Second DNA sequence")
    args = parser.parse_args()

    try:
        seq1 = args.seq1.strip().upper()
        seq2 = args.seq2.strip().upper()

        if not seq1 or not seq2:
            raise ValueError("Sequences must be non-empty strings.")

        distance = calculate_hamming_distance(seq1, seq2)
        print(distance)
        return 0
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
