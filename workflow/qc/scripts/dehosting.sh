#!/usr/bin/env bash
set -euo pipefail

##############################################
# Default parameters
##############################################
MAPQ=30
INC_FLAG=1
EXC_FLAG=12
MINLEN=50
MINFRAC=0.5

##############################################
# Usage
##############################################
usage() {
  cat <<EOF
Usage: $0 -i input.sam -o output.ids [options]

Required:
  -i FILE   Input SAM/BAM file
  -o FILE   Output read ID file

Optional:
  -q INT    Minimum MAPQ (default: 30)
  -f INT    Include SAM flag (default: 1 = read is paired)
  -F INT    Exclude SAM flag (default: 12 = -F 4 -F 8 -> Both mates are mapped, but we do not require “proper pair”)
  -l INT    Minimum aligned length in bp (default: 50)
  -r FLOAT  Minimum aligned fraction of read (default: 0.5)

Example:
  $0 -i sample.vs.Human.sam -o human_contam.ids
  $0 -i sample.bam -o human_contam.ids -q 40 -l 75 -r 0.7 -f 2
EOF
  exit 1
}

##############################################
# Parse arguments
##############################################
while getopts "i:o:q:f:F:l:r:h" opt; do
  case "$opt" in
    i) INPUT="$OPTARG" ;;
    o) OUTPUT="$OPTARG" ;;
    q) MAPQ="$OPTARG" ;;
    f) REQ_FLAG="$OPTARG" ;;
    F) EXC_FLAG="$OPTARG" ;;
    l) MINLEN="$OPTARG" ;;
    r) MINFRAC="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done

[[ -z "${INPUT:-}" || -z "${OUTPUT:-}" ]] && usage

##############################################
# Main logic
##############################################
samtools view -q "$MAPQ" -f "$INC_FLAG" -F "$EXC_FLAG" "$INPUT" |
awk -v minlen="$MINLEN" -v minfrac="$MINFRAC" '
BEGIN { FS = "\t" }

{
  # Read ID (same for both mates)
  read = $1

  # Total read length (from sequence field)
  seq_len = length($10)

  # Extract the CIGAR string
  cigar = $6

  # Variable to store total number of aligned bases (M)
  m = 0

  # Loop over every "number + M" pattern in the CIGAR string
  # Example: 63S21M -> finds "21M"
  while (match(cigar, /([0-9]+)M/)) {

    # Extract the numeric part (drop the 'M') and add to total
    m += substr(cigar, RSTART, RLENGTH - 1)

    # Remove the processed part of the CIGAR string
    # so we can find additional M blocks if they exist
    cigar = substr(cigar, RSTART + RLENGTH)
  }

  # Keep this read only if:
  # 1) aligned length ≥ minimum threshold
  # 2) aligned fraction ≥ minimum threshold
  if (m >= minlen && (m / seq_len) >= minfrac)
    good[read]++
}

END {
  for (r in good)
    if (good[r] == 2)
      print r
}
' > "$OUTPUT"