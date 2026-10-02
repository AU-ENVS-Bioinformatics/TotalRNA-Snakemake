#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  bash extract_glacial.sh ORDER_TAXIDS.txt [SAMPLE_DIR]

Arguments:
  ORDER_TAXIDS.txt  Overlapping SILVA order taxids, one per line.
  SAMPLE_DIR       Directory containing LARS sample folders (default: .).

Expected files within each sample folder:
  SAMPLE.silva.kraken
  SAMPLE.silva.report
  SAMPLE.prokaryotes_R1.fastq
  SAMPLE.prokaryotes_R2.fastq

Outputs within each sample folder:
  SAMPLE.glacial_orders_R1.fastq
  SAMPLE.glacial_orders_R2.fastq

Example:
  bash extract_glacial.sh glacier_silva_order_taxids.txt
EOF
}

if [[ $# -lt 1 || $# -gt 2 ]]; then
    usage >&2
    exit 1
fi

taxids_file="$1"
sample_dir="${2:-.}"
extractor=/data_2/rasmus/software/KrakenTools/extract_kraken_reads.py

[[ -r "$extractor" ]] || {
    echo "Cannot read extractor: $extractor" >&2
    exit 1
}

[[ -r "$taxids_file" && -s "$taxids_file" ]] || {
    echo "Taxid file is unreadable or empty: $taxids_file" >&2
    exit 1
}

[[ -d "$sample_dir" ]] || {
    echo "Cannot find sample directory: $sample_dir" >&2
    exit 1
}

# Validate taxids before processing any samples.
if ! awk '
    {
        sub(/\r$/, "")
        if ($0 !~ /^[0-9]+$/) {
            printf "Invalid taxid at line %d: %s\n", NR, $0 > "/dev/stderr"
            invalid = 1
        }
    }
    END { exit invalid }
' "$taxids_file"; then
    exit 1
fi

mapfile -t taxids < <(
    sed 's/\r$//' "$taxids_file" | LC_ALL=C sort -u
)

printf 'Extracting reads belonging to %s overlapping orders.\n' \
    "${#taxids[@]}" >&2

sample_count=0
shopt -s nullglob

for sample_path in "$sample_dir"/LARS[0-9]*; do
    [[ -d "$sample_path" ]] || continue

    sample="${sample_path##*/}"
    [[ "$sample" =~ ^LARS[0-9]+$ ]] || continue

    kraken="$sample_path/$sample.silva.kraken"
    report="$sample_path/$sample.silva.report"
    read1="$sample_path/$sample.prokaryotes_R1.fastq"
    read2="$sample_path/$sample.prokaryotes_R2.fastq"

    for file in "$kraken" "$report" "$read1" "$read2"; do
        [[ -r "$file" ]] || {
            echo "Cannot read: $file" >&2
            exit 1
        }
    done

    echo "Extracting glacier-associated order reads from $sample" >&2

    python "$extractor" \
        -k "$kraken" \
        -r "$report" \
        -s1 "$read1" \
        -s2 "$read2" \
        -t "${taxids[@]}" \
        --include-children \
        --fastq-output \
        --max 10000000000 \
        -o "$sample_path/$sample.glacial_orders_R1.fastq" \
        -o2 "$sample_path/$sample.glacial_orders_R2.fastq"

    sample_count=$((sample_count + 1))
done

if [[ "$sample_count" -eq 0 ]]; then
    echo "No LARS sample folders found in: $sample_dir" >&2
    exit 1
fi

printf 'Completed extraction for %s samples.\n' "$sample_count" >&2