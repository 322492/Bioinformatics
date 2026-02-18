#!/usr/bin/env python3
import sys


def reverse_complement(seq: str) -> str:
    complement = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
        "N": "N",
    }

    seq = seq.upper()
    rev = seq[::-1]

    out = []
    for base in rev:
        if base not in complement:
            raise ValueError(f"Invalid nucleotide: {base}")
        out.append(complement[base])

    return "".join(out)


def read_fasta(filepath: str):
    header = None
    sequence_parts = []

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line
            else:
                sequence_parts.append(line)

    sequence = "".join(sequence_parts)
    return header, sequence


def main():
    if len(sys.argv) != 2:
        print("Usage: python revcomp.py <input_fasta>", file=sys.stderr)
        sys.exit(2)

    infile = sys.argv[1]
    header, seq = read_fasta(infile)
    if header is None:
        print("Error: no FASTA header found.", file=sys.stderr)
        sys.exit(1)

    rc = reverse_complement(seq)

    print(header)
    print(rc)


if __name__ == "__main__":
    main()
