#!/usr/bin/env python3
"""
Exercise 3: sequence distance metrics.

Calculates:
1) Hamming distance
2) Composition vectors [A, T, G, C]
3) Manhattan distance on composition vectors
4) Euclidean distance on composition vectors
"""

from __future__ import annotations

import math
from collections import Counter
from typing import Dict, List, Sequence, Tuple


NUCLEOTIDE_ORDER: Tuple[str, ...] = ("A", "T", "G", "C")


def validate_sequence(seq: str) -> str:
    """Normalize and validate a DNA sequence."""
    normalized = seq.strip().upper()
    invalid = sorted(set(normalized) - set(NUCLEOTIDE_ORDER))
    if invalid:
        raise ValueError(
            f"Invalid symbols in sequence '{seq}': {', '.join(invalid)}. "
            "Allowed symbols: A, T, G, C."
        )
    return normalized


def hamming_distance(seq1: str, seq2: str) -> Tuple[int, List[Tuple[int, str, str]]]:
    """Return Hamming distance and list of mismatches as (position, s1, s2)."""
    if len(seq1) != len(seq2):
        raise ValueError(
            "Hamming distance requires sequences of equal length "
            f"(got {len(seq1)} and {len(seq2)})."
        )

    mismatches: List[Tuple[int, str, str]] = []
    for idx, (a, b) in enumerate(zip(seq1, seq2), start=1):
        if a != b:
            mismatches.append((idx, a, b))
    return len(mismatches), mismatches


def composition_vector(seq: str, order: Sequence[str] = NUCLEOTIDE_ORDER) -> Dict[str, int]:
    """Return nucleotide counts in a fixed order."""
    counts = Counter(seq)
    return {nt: counts.get(nt, 0) for nt in order}


def manhattan_distance(v1: Dict[str, int], v2: Dict[str, int], order: Sequence[str]) -> int:
    """Compute L1 distance between composition vectors."""
    return sum(abs(v1[nt] - v2[nt]) for nt in order)


def euclidean_distance(v1: Dict[str, int], v2: Dict[str, int], order: Sequence[str]) -> float:
    """Compute L2 distance between composition vectors."""
    return math.sqrt(sum((v1[nt] - v2[nt]) ** 2 for nt in order))


def format_comparison_rows(seq1: str, seq2: str) -> List[str]:
    """Build aligned rows for position-by-position sequence comparison."""
    rows: List[str] = []
    header = f"{'Pos':>3} | {'Seq1':^4} | {'Seq2':^4} | {'Match':^7}"
    sep = "-" * len(header)
    rows.append(header)
    rows.append(sep)
    for i, (a, b) in enumerate(zip(seq1, seq2), start=1):
        mark = "Yes" if a == b else "No"
        rows.append(f"{i:>3} | {a:^4} | {b:^4} | {mark:^7}")
    return rows


def format_vector(vec: Dict[str, int], order: Sequence[str]) -> str:
    """Format vector as [A:x, T:y, G:z, C:w]."""
    return "[" + ", ".join(f"{nt}:{vec[nt]}" for nt in order) + "]"


def main() -> None:
    # Set the two sequences here (exercise input or your custom test case).
    seq1_raw = "ATGCTAG"
    seq2_raw = "ATCCTAG"

    seq1 = validate_sequence(seq1_raw)
    seq2 = validate_sequence(seq2_raw)

    ham_dist, mismatches = hamming_distance(seq1, seq2)
    vec1 = composition_vector(seq1, NUCLEOTIDE_ORDER)
    vec2 = composition_vector(seq2, NUCLEOTIDE_ORDER)
    man_dist = manhattan_distance(vec1, vec2, NUCLEOTIDE_ORDER)
    euc_dist = euclidean_distance(vec1, vec2, NUCLEOTIDE_ORDER)

    print("=" * 64)
    print("Exercise 3 - Sequence Distance Metrics")
    print("=" * 64)
    print(f"Seq1: {seq1} (length={len(seq1)})")
    print(f"Seq2: {seq2} (length={len(seq2)})")
    print()

    print("Position-by-position comparison")
    for line in format_comparison_rows(seq1, seq2):
        print(line)
    print()

    print("3.1 Hamming Distance")
    print(f"Result: {ham_dist}")
    if mismatches:
        mismatch_text = ", ".join(
            f"pos {pos}: {a}->{b}" for pos, a, b in mismatches
        )
        print(f"Mismatches: {mismatch_text}")
    else:
        print("Mismatches: none")
    print()

    print("3.2 Composition Vectors [A, T, G, C]")
    print(f"Seq1 vector: {format_vector(vec1, NUCLEOTIDE_ORDER)}")
    print(f"Seq2 vector: {format_vector(vec2, NUCLEOTIDE_ORDER)}")
    print()

    print("3.3 Manhattan Distance (L1)")
    print(f"Result: {man_dist}")
    print()

    print("3.4 Euclidean Distance (L2)")
    print(f"Result: {euc_dist:.6f}")
    print()


if __name__ == "__main__":
    main()
