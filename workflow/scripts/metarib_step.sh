#!/bin/bash

##############################################################################
# metarib_step.sh
# 
# Standalone script to run a single MetaRib iteration cycle:
# 1. EMIRGE amplicon assembly on subsampled reads
# 2. Dereplication of assembled contigs
# 3. BBMap alignment to find unmapped reads for next iteration
#
# Usage: metarib_step.sh [options]
##############################################################################

set -euo pipefail

# Default values
N=10000
THREADS=8
MAP_PARA="minid=0.96 maxindel=1 minhits=2 idfilter=0.98"
CLS_PARA="fo=t ow=t c=t mcs=1 e=5 mid=99"
SAMPLE="sample"
ITER=0

# Function to display usage
usage() {
    cat << EOF
Usage: $0 -1 <r1_file1> <r1_file2> ... -2 <r2_file1> <r2_file2> ... [options]

Required arguments:
  -1, --r1              One or more R1 reads fastq files (space-separated)
  -2, --r2              One or more R2 reads fastq files (space-separated)
  -c, --contigs         Path to input contigs fasta file
  -o, --output-dir      Output directory for results
  -f, --ref-db          Reference database for EMIRGE (FASTA file)
  -b, --bt-idx          Bowtie2 index prefix for reference database
  -e, --em-para         EMIRGE parameters


Optional arguments:
  -s, --sample          Sample name (default: 'sample')
  -i, --iter            Iteration number (default: 0)
  -n, --num-reads       Number of reads to subsample (default: 10000)
  -t, --threads         Number of threads (default: 8)
  --map-para            BBMap parameters (default: '$MAP_PARA')
  --cls-para            Dedupe parameters (default: '$CLS_PARA')
  -h, --help            Display this help message

EOF
    exit 1
}

# Parse arguments - collect multiple files for -1 and -2
R1=()
R2=()
CONTIGS=""
OUTPUT_DIR=""
REF_DB=""
BT_IDX=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -1|--r1) 
            shift
            while [[ $# -gt 0 && "$1" != -* ]]; do
                R1+=("$1")
                shift
            done
            ;;
        -2|--r2) 
            shift
            while [[ $# -gt 0 && "$1" != -* ]]; do
                R2+=("$1")
                shift
            done
            ;;
        -c|--contigs) CONTIGS="$2"; shift 2 ;;
        -o|--output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        -f|--ref-db) REF_DB="$2"; shift 2 ;;
        -b|--bt-idx) BT_IDX="$2"; shift 2 ;;
        -s|--sample) SAMPLE="$2"; shift 2 ;;
        -i|--iter) ITER="$2"; shift 2 ;;
        -n|--num-reads) N="$2"; shift 2 ;;
        -t|--threads) THREADS="$2"; shift 2 ;;
        --em-para) EM_PARA="$2"; shift 2 ;;
        --map-para) MAP_PARA="$2"; shift 2 ;;
        --cls-para) CLS_PARA="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

# Validate required arguments
if [[ ${#R1[@]} -eq 0 || ${#R2[@]} -eq 0 || -z "$CONTIGS" || -z "$OUTPUT_DIR" || -z "$REF_DB" || -z "$BT_IDX" ]]; then
    echo "Error: Missing required arguments" >&2
    usage
fi

# Check R1 and R2 have same number of files
if [[ ${#R1[@]} -ne ${#R2[@]} ]]; then
    echo "Error: Number of R1 and R2 files must be equal" >&2
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

dirname=$(realpath "$(dirname "${R1[0]}")")

# Convert relative paths to absolute paths before changing directory
for i in "${!R1[@]}"; do
    R1[$i]=$(realpath "${R1[$i]}")
    R2[$i]=$(realpath "${R2[$i]}")
done

CONTIGS=$(realpath "$CONTIGS")
REF_DB=$(realpath "$REF_DB")
# BT_IDX is a prefix, so we need to handle it specially
BT_IDX_DIR=$(dirname "$BT_IDX")
BT_IDX_BASE=$(basename "$BT_IDX")
BT_IDX="$(realpath "$BT_IDX_DIR")/$BT_IDX_BASE"
OUTPUT_DIR=$(realpath "$OUTPUT_DIR")
# Output next
CONTIGS_NEXT="$OUTPUT_DIR/contigs_derep_next.fasta"
R1_NEXT="$OUTPUT_DIR/unmapped_R1_next.fastq"
R2_NEXT="$OUTPUT_DIR/unmapped_R2_next.fastq"

# Parsing EMIRGE parameters into individual variables for clarity
# EM_PARA format: "--phred33 -l 125 -i 250 -s 50 -a 20 -n 20"
# Extract individual parameter values using grep and sed
EM_PARA_l=$(echo "$EM_PARA" | grep -oP '\-l\s+\K\d+' || echo "125")        # max_read_length
EM_PARA_i=$(echo "$EM_PARA" | grep -oP '\-i\s+\K\d+' || echo "250")        # insert_mean
EM_PARA_s=$(echo "$EM_PARA" | grep -oP '\-s\s+\K\d+' || echo "50")         # insert_stddev
EM_PARA_a=$(echo "$EM_PARA" | grep -oP '\-a\s+\K\d+' || echo "20")         # EMIRGE parameter a
EM_PARA_n=$(echo "$EM_PARA" | grep -oP '\-n\s+\K\d+' || echo "20")         # EMIRGE parameter n
EM_PARA_phred=$(echo "$EM_PARA" | grep -oP '\-\-phred\d+' || echo "--phred33")  # phred quality

##############################################################################
# MAIN EXECUTION
##############################################################################

echo "Started at: $(date)"
echo "Working directory: $OUTPUT_DIR"

# Step 1: EMIRGE amplicon assembly
echo ""
echo "Step 1: Running EMIRGE on $N subsampled reads..."

cd "$OUTPUT_DIR"

seed=$(( (RANDOM << 15 | RANDOM) % 999999 + 1 ))

cat "$dirname"/*.1.fq | seqtk sample -s"$seed" - "$N" > subsample_R1.fastq 2>&1
cat "$dirname"/*.2.fq | seqtk sample -s"$seed" - "$N" > subsample_R2.fastq 2>&1


max_read_length=$(awk 'NR%4==2 {print length($0)}' subsample_R1.fastq | sort -n | tail -1)

vsearch \
  --fastq_mergepairs subsample_R1.fastq \
  --reverse subsample_R2.fastq \
  --fastqout merged.fastq \
  --fastq_minovlen 20

awk 'NR%4==2 {print length($0)}' merged.fastq > merged_lengths.txt

mean_dist=$(awk '{sum+=$1} END {printf "%.0f\n", sum/NR}' merged_lengths.txt)

stddev_dist=$(awk '{
  sum+=$1; 
  sumsq+=$1*$1
} END {
  mean=sum/NR;
  stddev=sqrt(sumsq/NR - mean*mean);
  printf "%.0f\n", stddev
}' merged_lengths.txt)

emirge_amplicon.py emirge_subset \
    -1 subsample_R1.fastq -2 subsample_R2.fastq \
    --max_read_length "$max_read_length" --insert_mean "$mean_dist" --insert_stddev "$stddev_dist" --processors "$EM_PARA_a" --iterations "$EM_PARA_n" \
    "$EM_PARA_phred" --fasta_db "$REF_DB" --bowtie_db "$BT_IDX" 2>&1


# Step 2: Dereplication
echo ""
echo "Step 2: Deduplicating contigs..."

cat emirge_subset/iter.*/*.fasta "$CONTIGS" > contigs.combined.fasta

# sortbyname.sh in=contigs.combined.fasta out=contigs.sorted.fasta length descending 2>&1
# reformat.sh in=contigs.sorted.fasta out=contigs.formatted.fasta uniquenames 2>&1
# dedupe.sh in=contigs.formatted.fasta out="$CONTIGS_NEXT" outd=contigs.duplicates.fasta $CLS_PARA 2>&1

# Deduplicate and keep size info
vsearch --fasta_width 0 \
  --derep_fulllength contigs.combined.fasta \
  --output temp_derep.fasta \
  --relabel contig_

vsearch --sortbysize temp_derep.fasta \
  --output "$CONTIGS_NEXT"

contig_count=$(grep -c '^>' "$CONTIGS_NEXT" || echo 0)
echo "$contig_count" > iteration_report.txt
echo "Assembled contigs: $contig_count"

# Map reads and extract unmapped in parallel
echo ""
echo "Step 3: Mapping reads to assembled contigs..."
bbmap.sh ref="$CONTIGS_NEXT" 2>&1

mkdir -p "$OUTPUT_DIR"/unmapped_data

# Handle both single files and multiple input reads
r1_files=($R1)
r2_files=($R2)
num_samples=${#r1_files[@]}
threads_per_sample=$(( THREADS / num_samples ))
[ "$threads_per_sample" -lt 1 ] && threads_per_sample=1

for i in "${!r1_files[@]}"; do
    r1="${r1_files[$i]}"
    r2="${r2_files[$i]}"
    sample=$(basename "$r1" .1.fq)
    bbmap.sh in="$r1" in2="$r2" \
        outu=unmapped_data/${sample}.1.fq \
        outu2=unmapped_data/${sample}.2.fq \
        ref="$CONTIGS_NEXT" \
        # threads=2 2>&1 &
        threads=$threads_per_sample 2>&1 &
done


# Sum lines from all files
total_lines=$(wc -l unmapped_data/*.fq | tail -1 | awk '{print $1}')
total_reads=$(( total_lines / 4 / 2 ))  # Divide by 4 (FASTQ) and 2 (paired files)
echo "$total_reads" >> iteration_report.txt
echo "Unmapped reads remaining: $total_reads"


echo "Iteration complete at: $(date)"
echo "Summary: $contig_count contigs assembled"
echo "========================================"

# Create output summary
echo ""
echo "=========================================="
echo "Step completed successfully!"
echo "Output directory: $OUTPUT_DIR"
echo "Contigs: $CONTIGS_NEXT"
echo "Unmapped R1: $R1_NEXT"
echo "Unmapped R2: $R2"
echo "Duplicates: contigs.duplicates.fasta"
echo "Report: iteration_report.txt"
echo "========================================"

exit 0