#!/usr/bin/env bash
# Exercise 4.2 - Global alignment with EMBOSS needle
# This script:
# 1) generates seqA.fasta and seqB.fasta
# 2) runs needle with default scoring parameters
# 3) saves result to global_alignment.txt

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

if ! command -v needle >/dev/null 2>&1; then
  echo "Error: EMBOSS needle is not available in PATH."
  echo "Install EMBOSS first, for example:"
  echo "  conda install -n bioinfo_lab3 -c conda-forge emboss -y"
  exit 1
fi

cat > seqA.fasta <<'EOF'
>seqA
GATCTA
EOF

cat > seqB.fasta <<'EOF'
>seqB
GAACGTA
EOF

# -auto disables interactive prompts, keeping the run fully scripted.
needle \
  -asequence seqA.fasta \
  -bsequence seqB.fasta \
  -outfile global_alignment.txt \
  -auto

echo "Done. Generated files:"
echo "  - ${SCRIPT_DIR}/seqA.fasta"
echo "  - ${SCRIPT_DIR}/seqB.fasta"
echo "  - ${SCRIPT_DIR}/global_alignment.txt"
