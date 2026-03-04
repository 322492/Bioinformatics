#!/usr/bin/env python3
"""
Exercise 6: Generate all k-mers from a sequence.
"""

from __future__ import annotations


def get_kmers(sequence: str, k: int) -> list[str]:
    """
    Return all overlapping k-mers from the input sequence.

    Example:
    sequence = "GATCGATC", k = 3
    result   = ["GAT", "ATC", "TCG", "CGA", "GAT", "ATC"]
    """
    seq = sequence.strip().upper()

    if not seq:
        raise ValueError("Sequence must be a non-empty string.")
    if k <= 0:
        raise ValueError("k must be a positive integer.")
    if k > len(seq):
        raise ValueError(
            f"k cannot be larger than sequence length (k={k}, length={len(seq)})."
        )

    return [seq[i : i + k] for i in range(len(seq) - k + 1)]


def main() -> None:
    # Set input parameters here.
    sequence = "GATCGATC"
    k = 3

    print("=" * 60)
    print("Exercise 6 - k-mer generation")
    print("=" * 60)
    print(f"Input sequence : {sequence}")
    print(f"Window size (k): {k}")
    print()

    try:
        kmers = get_kmers(sequence, k)
        print("Output k-mers:")
        print(kmers)
        print()
        print(f"Total k-mers: {len(kmers)}")
    except ValueError as exc:
        print(f"Error: {exc}")


if __name__ == "__main__":
    main()
