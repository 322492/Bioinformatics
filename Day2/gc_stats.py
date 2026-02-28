#!/usr/bin/env python3
import sys

def emit_counts(a: int, t: int, g: int, c: int) -> None:
    length = a + t + g + c
    gc_pct = (100.0 * (g + c) / length) if length > 0 else 0.0
    # TSV: #A #T #G #C GC%
    sys.stdout.write(f"{a}\t{t}\t{g}\t{c}\t{gc_pct:.6f}\n")

def main() -> None:
    # Print header
    sys.stdout.write("#A\t#T\t#G\t#C\tGC%\n")

    a = t = g = c = 0
    in_seq = False

    total_g = total_c = total_len = 0

    for line in sys.stdin:
        if not line:
            continue
        line = line.strip()
        if not line:
            continue

        if line.startswith(">"):
            # flush previous sequence (if any)
            if in_seq:
                emit_counts(a, t, g, c)
                total_g += g
                total_c += c
                total_len += (a + t + g + c)
                a = t = g = c = 0
            in_seq = True
            continue

        # sequence line
        in_seq = True
        for ch in line:
            # FASTA could be lowercase; we normalize
            if ch == "A" or ch == "a":
                a += 1
            elif ch == "T" or ch == "t":
                t += 1
            elif ch == "G" or ch == "g":
                g += 1
            elif ch == "C" or ch == "c":
                c += 1
            else:
                # ignore ambiguous bases (N, etc.)
                pass

    # flush last sequence
    if in_seq:
        emit_counts(a, t, g, c)
        total_g += g
        total_c += c
        total_len += (a + t + g + c)

    weighted_gc = (100.0 * (total_g + total_c) / total_len) if total_len > 0 else 0.0
    sys.stdout.write(f"---\t{weighted_gc:.6f}\n")

if __name__ == "__main__":
    main()