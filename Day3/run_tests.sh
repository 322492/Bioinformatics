#!/usr/bin/env bash
# Exercise 5.3 - Automated tests for hamming_distance.py
# Generates hamming.csv from scratch.

set -u
set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

OUT_CSV="hamming.csv"
TMP_ERR="$(mktemp)"
trap 'rm -f "${TMP_ERR}"' EXIT

printf "pair,distance\n" > "${OUT_CSV}"

run_pair() {
  local pair_name="$1"
  local seq1="$2"
  local seq2="$3"
  local distance

  if distance="$(python3 hamming_distance.py --seq1 "${seq1}" --seq2 "${seq2}" 2>"${TMP_ERR}")"; then
    printf "%s,%s\n" "${pair_name}" "${distance}" >> "${OUT_CSV}"
  else
    # Required behavior for failing input (e.g., Pair C): record -1.
    printf "%s,-1\n" "${pair_name}" >> "${OUT_CSV}"
    printf "Warning: %s failed (%s)\n" "${pair_name}" "$(cat "${TMP_ERR}")" >&2
  fi
}

run_pair "Pair A" "ATGCGTACGTAGCTA" "ATGCCTACGTAGCTA"
run_pair "Pair B" "GAGCCTACTAACGGGAT" "CATCGTAATGACGGCCT"
run_pair "Pair C" "ATGC" "ATG"

echo "Done. Generated ${SCRIPT_DIR}/${OUT_CSV}"
