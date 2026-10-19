"""Expected completion outputs for viral sequence analysis."""

from pathlib import Path
from typing import List, Sequence


def viral_output(
    results_dir: Path,
    samples: Sequence[str],
    config: dict,
) -> List[str]:
    """Return files required to complete viral sequence analysis."""
    viral_dir = results_dir / "Viral_Sequence_Analysis"

    genomad_dir = (
        viral_dir
        / "identification"
        / "genomad"
        / "transcripts.filtered_summary"
    )

    outputs: List[str] = [
        f"{viral_dir}/transcript/transcripts.filtered.stats.tsv",
        f"{genomad_dir}/transcripts.filtered_plasmid_summary.tsv",
        f"{genomad_dir}/transcripts.filtered_virus_summary.tsv",
        f"{genomad_dir}/transcripts.filtered_virus.fna",
        f"{genomad_dir}/transcripts.filtered_virus_genes.tsv",
        f"{viral_dir}/quality/checkv/quality_summary.tsv",
    ]

    if config["viral_sequence_analysis"].get("postprocessing", True):
        outputs += [
            f"{viral_dir}/final_result/viral_candidates.tsv",
            f"{viral_dir}/final_result/viral_eggnog_annotations.tsv",
        ]

        outputs += [
            f"{viral_dir}/abundance/{sample}/viral_quant.sf"
            for sample in samples
        ]

    return outputs