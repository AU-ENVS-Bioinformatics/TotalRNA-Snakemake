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
Usage: $0 -1 <r1_file> -2 <r2_file> -c <contigs_file> -o <output_dir> -f <ref_db> -b <bt_idx> [options]

Required arguments:
  -1, --r1              Path to forward (R1) reads fastq file
  -2, --r2              Path to reverse (R2) reads fastq file
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

# Parse arguments
R1=""
R2=""
CONTIGS=""
OUTPUT_DIR=""
REF_DB=""
BT_IDX=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -1|--r1) R1="$2"; shift 2 ;;
        -2|--r2) R2="$2"; shift 2 ;;
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
if [[ -z "$R1" || -z "$R2" || -z "$CONTIGS" || -z "$OUTPUT_DIR" || -z "$REF_DB" || -z "$BT_IDX" ]]; then
    echo "Error: Missing required arguments" >&2
    usage
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Convert relative paths to absolute paths before changing directory
R1=$(realpath "$R1")
R2=$(realpath "$R2")
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

# Verify input files exist
for file in "$R1" "$R2" "$CONTIGS" "$REF_DB"; do
    if [[ ! -f "$file" ]]; then
        echo "Error: Input file not found: $file" >&2
        exit 1
    fi
done

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
seqtk sample -s"$seed" "$R1" "$N" > subsample_R1.fastq 2>&1
seqtk sample -s"$seed" "$R2" "$N" > subsample_R2.fastq 2>&1


# awk 'NR%4==2 {print length($0)}' results/metarib/work/iter.1/subsample_R1.fastq | sort -n | tail


emirge_amplicon.py emirge_subset \
    -1 subsample_R1.fastq -2 subsample_R2.fastq \
    --max_read_length "$EM_PARA_l" --insert_mean "$EM_PARA_i" --insert_stddev "$EM_PARA_s" --processors "$EM_PARA_a" --iterations "$EM_PARA_n" \
    "$EM_PARA_phred" --fasta_db "$REF_DB" --bowtie_db "$BT_IDX" \
    2>&1

# emirge_amplicon.py emirge_subset \
#     -1 subsample_R1.fastq -2 subsample_R2.fastq \
#     --phred33 -l 125 -i 250 -s 50 -a 20 -n 10 \
#     --fasta_db "$REF_DB" --bowtie_db "$BT_IDX" 2>&1

 
# Step 2: Dereplication
echo ""
echo "Step 2: Deduplicating contigs..."

cat emirge_subset/iter.*/*.fasta "$CONTIGS" > contigs.combined 2>&1

sortbyname.sh in=contigs.combined out=contigs.sorted length descending 2>&1
reformat.sh in=contigs.sorted out=contigs.formatted uniquenames 2>&1
dedupe.sh in=contigs.formatted out="$CONTIGS_NEXT" outd=contigs.duplicates.fasta $CLS_PARA 2>&1

contig_count=$(grep -c '^>' "$CONTIGS_NEXT" || echo 0)
echo "$contig_count" > iteration_report.txt
echo "Assembled contigs: $contig_count"

# Step 3: Map reads and extract unmapped
echo ""
echo "Step 3: Mapping reads to assembled contigs..."

bbmap.sh ref="$CONTIGS_NEXT" in1="$R1" in2="$R2" \
    threads=$THREADS $MAP_PARA outu=unmapped.fq ow=t \
    statsfile=bbmap.stats.txt sortscafs=t \
    scafstats=bbmap.scafstats.txt covstats=bbmap.covstats.txt 2>&1

reformat.sh in=unmapped.fq out1="$R1_NEXT" out2="$R2_NEXT" 2>&1

r1_reads=$(( $(wc -l < "$R1_NEXT") / 4 ))
echo "$r1_reads" >> iteration_report.txt
echo "Unmapped reads remaining: $r1_reads"

echo ""
echo "Iteration complete at: $(date)"
echo "Summary: $contig_count contigs assembled, $r1_reads reads remain unmapped"
echo "========================================"

# Create output summary
echo ""
echo "=========================================="
echo "Step completed successfully!"
echo "Output directory: $OUTPUT_DIR"
echo "Contigs: $CONTIGS"
echo "Unmapped R1: $R1"
echo "Unmapped R2: $R2"
echo "Duplicates: contigs.duplicates.fasta"
echo "Report: iteration_report.txt"
echo "=========================================="

exit 0
