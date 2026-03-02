import sys
import pysam
import pandas as pd
from collections import defaultdict
from typing import Dict
from lib.utils import regex_filename

files = snakemake.input.bam
output = snakemake.output.tsv
log = snakemake.log

def contig_mapped_read_length(bam_path: str) -> pd.DataFrame:
    """Compute total mapped read length per contig from a BAM file."""
    contig_totals: Dict[str, int] = defaultdict(int)
    print(f"Processing BAM file: {bam_path}", file=sys.stderr)
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            contig = bam.get_reference_name(read.reference_id)
            contig_totals[contig] += read.query_alignment_length

    df = pd.DataFrame(
        [(contig, total) for contig, total in contig_totals.items()],
        columns=["contig", "read_length"],
    )
    print(f"Finished processing BAM file: {bam_path}", file=sys.stderr)
    df["sample"] = regex_filename(bam_path)

    return df[["sample", "contig", "read_length"]]

with open(log[0], "w") as f:
    sys.stderr = sys.stdout = f
    print("Calculating total mapped read lengths per contig...", file=sys.stderr)
    print(files, file=sys.stderr)
    df = contig_mapped_read_length(files)
    print("Finished calculating total mapped read lengths per contig.", file=sys.stderr)
    print(f"Writing output to {output}...", file=sys.stderr)
    df.to_csv(output, sep="\t", index=False)
    print(f"Finished writing output to {output}.", file=sys.stderr)