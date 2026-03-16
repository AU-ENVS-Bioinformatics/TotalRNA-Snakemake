
set -euo pipefail



# Function to display usage
usage() {
    cat << EOF
Usage: $0 --dirname <input_directory> --num_reads <number_of_reads>
Required arguments:
  --dirname            Directory containing input .1.fq and .2.fq files
  --num_reads          Number of reads to subsample
EOF
    exit 1
}


while [[ $# -gt 0 ]]; do
    case $1 in
        --dirname|-d) dirname="$2"; shift 2 ;;
        --num_reads|-n ) N="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done


# Validate required arguments
if [[ -z "${dirname:-}" || -z "${N:-}" ]]; then
    echo "Error: Missing required arguments" >&2
    usage
fi

pushd "$dirname"

seed=$(( (RANDOM << 15 | RANDOM) % 999999 + 1 ))

cat *.1.fq | seqtk sample -s"$seed" - "$N" > subsample_R1.fastq 2>&1
cat *.2.fq | seqtk sample -s"$seed" - "$N" > subsample_R2.fastq 2>&1


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

rm merged.fastq merged_lengths.txt

popd


EM_PARA="--phred33 -l $max_read_length -i $mean_dist -s $stddev_dist -a 40 -n 10"