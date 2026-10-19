from typing import List
from pathlib import Path
from snakemake.io import expand
from utils.reconstruction import validate_reconstruction_config

def reconstruction_outputs(
    results_dir: Path,
    samples: List[str],
    config: dict,
    active_kingdoms: List[str],
) -> List[str]:
    """
    Return all expected reconstruction output files.
    """

    validate_reconstruction_config(config)

    secondary_method = config["RNA_classification"].get(
        "secondary_method", "none"
    )

    outputs: List[str] = []

    # Sample-specific SSU reconstruction
    outputs += expand(
        f"{results_dir}/Assembly/{{sample}}/SSU/rnaspades/transcripts.fasta",
        sample=samples,
    )

    # Combined SSU reconstruction
    outputs.append(
        f"{results_dir}/Assembly/concatenated/SSU/cross_sample_SSU.fasta"
    )

    # rRNA predictions for each configured Barrnap kingdom
    outputs += expand(
        f"{results_dir}/Assembly/SSU_candidates/barrnap/{{kingdom}}/cross_sample_SSU_{{kingdom}}.barrnap.gff",
        kingdom=config["reconstruction"]["barrnap"]["kingdom"],
    )

    outputs.append(f"{results_dir}/Assembly/SSU_candidates/barrnap/barrnap_kingdom_status.tsv")

    outputs += expand(
        f"{results_dir}/Assembly/SSU_candidates/postprocessing/{{kingdom}}/SSU_{{kingdom}}_filtered.fasta",
        kingdom=active_kingdoms,
    )

    outputs += expand(
        f"{results_dir}/Assembly/SSU_candidates/postprocessing/{{kingdom}}/SSU_{{kingdom}}_dereplicated.fasta",
        kingdom=active_kingdoms,
    )

    outputs += expand(
        f"{results_dir}/Assembly/SSU_candidates/postprocessing/{{kingdom}}/SSU_{{kingdom}}_nonchimeras.fasta",
        kingdom=active_kingdoms,
    )

    outputs += expand(
        f"{results_dir}/Assembly/SSU_candidates/postprocessing/{{kingdom}}/SSU_{{kingdom}}_catalogue.fasta",
        kingdom=active_kingdoms,
    )
    
    # Optional sample-specific ITS-candidate reconstruction
    if secondary_method == "ITS":
        outputs += expand(
            f"{results_dir}/Assembly/{{sample}}/ITS_candidates/rnaspades/transcripts.fasta",
            sample=samples,
        )

        outputs.append(
            f"{results_dir}/Assembly/concatenated/ITS/cross_sample_ITS.fasta"
        )

    # Cross-sample non-rRNA reconstruction
    outputs.append(
        f"{results_dir}/Assembly/coassembly/non_rRNA/rnaspades/transcripts.fasta"
    )

    return outputs