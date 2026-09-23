from typing import List
from pathlib import Path

def viral_output(
    results_dir: Path,
    config: dict,
) -> List[str]:
    """
    Return expected viral sequence analysis outputs.
    """

    viral_dir = f"{results_dir}/Viral_Sequence_Analysis"

    outputs: List[str] = []

    # Standardized final ORF outputs
    outputs += [
        f"{viral_dir}/transcript/transcripts.filtered.fasta",
        f"{viral_dir}/identification/genomad/transcripts_summary/transcripts_plasmid_summary.tsv",
        f"{viral_dir}/identification/genomad/transcripts_summary/transcripts_virus_summary.tsv",
        f"{viral_dir}/quality/checkv/quality_summary.tsv",
    ]

    return outputs