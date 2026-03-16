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

emirge_amplicon.py emirge_subset \
    -1 subsample_R1.fastq -2 subsample_R2.fastq \
    --max_read_length "$EM_PARA_l" --insert_mean "$EM_PARA_i" --insert_stddev "$EM_PARA_s" --processors "$EM_PARA_a" --iterations "$EM_PARA_n" \
    "$EM_PARA_phred" --fasta_db "$REF_DB" --bowtie_db "$BT_IDX" 2>&1


# Step 2: Dereplication
echo ""
echo "Step 2: Deduplicating contigs..."
cat emirge_subset/iter.*/iter.*.cons.fasta "$CONTIGS" > contigs.combined.fasta

# Deduplicate info - keep original labels and add ;size= suffix
vsearch \
  --derep_fulllength contigs.combined.fasta \
  --output contigs_derep_next.fasta  \
  --sizeout 2>&1

sed -i 's/;size=/_/g' contigs_derep_next.fasta 

contig_count=$(grep -c '^>' contigs_derep_next.fasta || echo 0)
echo "$contig_count" > iteration_report.txt
echo "Assembled contigs: $contig_count"

# Map reads and extract unmapped in parallel
echo ""
echo "Step 3: Mapping reads to assembled contigs..."
bbmap.sh -Xmx5g ref=contigs_derep_next.fasta 2>&1

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
    bbmap.sh -Xmx5g in="$r1" in2="$r2" \
        outu=unmapped_data/${sample}.1.fq \
        outu2=unmapped_data/${sample}.2.fq \
        ref=contigs_derep_next.fasta \
        "$MAP_PARA" \
        ow=t \
        threads=$threads_per_sample 2>&1
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
echo "Contigs: contigs_derep_next.fasta"
echo "Report: iteration_report.txt"
echo "========================================"

exit 0