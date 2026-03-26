#!/usr/bin/env python3
from Bio import Align
from Bio.Align import substitution_matrices


def main():
    seq1 = "HEAGAWGHEE"
    seq2 = "PAWHEAE"
    matrix = substitution_matrices.load("BLOSUM62")

    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = matrix
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    aligner.mode = "global"
    global_score = aligner.score(seq1, seq2)

    aligner.mode = "local"
    local_score = aligner.score(seq1, seq2)

    with open("Day6/exercise4_results.txt", "w", encoding="utf-8") as handle:
        handle.write(f"global_score={global_score}\n")
        handle.write(f"local_score={local_score}\n")


if __name__ == "__main__":
    main()
