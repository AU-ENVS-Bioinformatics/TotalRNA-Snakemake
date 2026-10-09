from pathlib import Path
from typing import List, Sequence


def viral_output(
    results_dir: Path,
    samples: Sequence[str],
    config: dict
) -> List[str]:
    """Return files required for completion of the viral analysis module."""
    viral_config = config.get("viral_sequence_analysis", {})
    if not viral_config.get("enabled", True):
        return []

    viral_dir = f"{results_dir}/Viral_Sequence_Analysis"
    genomad_dir = f"{viral_dir}/identification/genomad/transcripts.filtered_summary"
    outputs: List[str] = [
        f"{viral_dir}/transcript/transcripts.filtered.stats.tsv",
        f"{genomad_dir}/transcripts.filtered_plasmid_summary.tsv",
        f"{genomad_dir}/transcripts.filtered_virus_summary.tsv",
        f"{genomad_dir}/transcripts.filtered_virus.fna",
        f"{genomad_dir}/transcripts.filtered_virus_genes.tsv",
        f"{viral_dir}/quality/checkv/quality_summary.tsv",
    ]
    if not viral_config.get("postprocessing", True):
        return outputs

    outputs += [
        f"{viral_dir}/final_result/viral_candidates.tsv",
        f"{viral_dir}/final_result/viral_eggnog_annotations.tsv",
    ]
    outputs += [
        f"{viral_dir}/abundance/{sample}/viral_quant.sf"
        for sample in samples
    ]
    return outputs
