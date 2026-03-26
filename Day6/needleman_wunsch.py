import argparse


def needleman_wunsch(seq1, seq2, match_score, mismatch_score, gap_score):
    rows = len(seq1) + 1
    cols = len(seq2) + 1

    dp = [[0 for _ in range(cols)] for _ in range(rows)]
    tb = [[None for _ in range(cols)] for _ in range(rows)]

    for i in range(1, rows):
        dp[i][0] = dp[i - 1][0] + gap_score
        tb[i][0] = "up"
    for j in range(1, cols):
        dp[0][j] = dp[0][j - 1] + gap_score
        tb[0][j] = "left"

    for i in range(1, rows):
        for j in range(1, cols):
            diag = dp[i - 1][j - 1] + (
                match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score
            )
            up = dp[i - 1][j] + gap_score
            left = dp[i][j - 1] + gap_score

            best = max(diag, up, left)
            dp[i][j] = best

            if best == diag:
                tb[i][j] = "diag"
            elif best == up:
                tb[i][j] = "up"
            else:
                tb[i][j] = "left"

    return dp, tb


def traceback_alignment(seq1, seq2, tb):
    i = len(seq1)
    j = len(seq2)

    aln1 = []
    aln2 = []
    path = [(i, j)]

    while i > 0 or j > 0:
        move = tb[i][j]
        if move == "diag":
            aln1.append(seq1[i - 1])
            aln2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif move == "up":
            aln1.append(seq1[i - 1])
            aln2.append("-")
            i -= 1
        elif move == "left":
            aln1.append("-")
            aln2.append(seq2[j - 1])
            j -= 1
        else:
            break
        path.append((i, j))

    aln1.reverse()
    aln2.reverse()
    return "".join(aln1), "".join(aln2), path


def format_dp_matrix(seq1, seq2, dp):
    col_labels = ["-"] + list(seq2)
    row_labels = ["-"] + list(seq1)

    width = max(
        3,
        max(len(str(cell)) for row in dp for cell in row),
        max(len(label) for label in col_labels + row_labels),
    )

    lines = []
    header = f"{'':>{width}} " + " ".join(f"{c:>{width}}" for c in col_labels)
    lines.append(header)

    for i, row in enumerate(dp):
        line = f"{row_labels[i]:>{width}} " + " ".join(f"{v:>{width}}" for v in row)
        lines.append(line)

    return "\n".join(lines)


def format_traceback_path(path):
    return " -> ".join(f"({i},{j})" for i, j in path)


def main():
    parser = argparse.ArgumentParser(
        description="Needleman-Wunsch global alignment with DP matrix and traceback output."
    )
    parser.add_argument("--seq1", required=True, help="First sequence")
    parser.add_argument("--seq2", required=True, help="Second sequence")
    parser.add_argument("--match", type=int, default=1, help="Match score (default: 1)")
    parser.add_argument(
        "--mismatch", type=int, default=-1, help="Mismatch score (default: -1)"
    )
    parser.add_argument("--gap", type=int, default=-1, help="Gap score (default: -1)")
    args = parser.parse_args()

    seq1 = args.seq1.upper()
    seq2 = args.seq2.upper()

    dp, tb = needleman_wunsch(seq1, seq2, args.match, args.mismatch, args.gap)
    aligned1, aligned2, path = traceback_alignment(seq1, seq2, tb)

    print("DP scoring matrix:")
    print(format_dp_matrix(seq1, seq2, dp))
    print()
    print("Traceback path (from bottom-right to top-left):")
    print(format_traceback_path(path))
    print()
    print("Final alignment score:")
    print(dp[len(seq1)][len(seq2)])
    print()
    print("Optimal alignment:")
    print(aligned1)
    print(aligned2)


if __name__ == "__main__":
    main()
