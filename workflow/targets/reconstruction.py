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
    Return reconstruction outputs, optionally including Barrnap predictions and SSU postprocessing.
    """

    validate_reconstruction_config(config)

    assembly_refinement = config["reconstruction"].get("refinement","coassembly")

    secondary_method = config["RNA_classification"].get(
        "secondary_method", "none"
    )

    ssu_prediction = config["reconstruction"]["ssu_postprocessing"].get(
        "prediction_method", True
    )

    if not isinstance(ssu_prediction, bool):
        raise ValueError(
            "reconstruction.ssu_postprocessing.prediction_method "
            "must be true or false, without quotation marks."
        )

    outputs: List[str] = []

    # Separate cross-sample co-assemblies for each configured read category.
    if assembly_refinement in ["coassembly", "both", "all"]:
        coassembly_targets = config["reconstruction"].get("coassembly_targets", ["non_rRNA"])

        outputs += expand(
            f"{results_dir}/Assembly/coassembly/{{target}}/rnaspades/transcripts.fasta",
            target=coassembly_targets,
        )

    if assembly_refinement in ["samples","both","all"]:
        # Sample-specific SSU reconstruction
        outputs += expand(
            f"{results_dir}/Assembly/{{sample}}/SSU/rnaspades/transcripts.fasta",
            sample=samples,
        )

        # Combined SSU reconstruction
        outputs.append(
            f"{results_dir}/Assembly/concatenated/SSU/cross_sample_SSU.fasta"
        )

    # Optional Barrnap predictions and SSU postprocessing.
    if ssu_prediction:
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
    if (assembly_refinement in ["samples", "both", "all"] and secondary_method == "ITS"):
        outputs += expand(
            f"{results_dir}/Assembly/{{sample}}/ITS_candidates/rnaspades/transcripts.fasta",
            sample=samples,
        )

        outputs.append(
            f"{results_dir}/Assembly/concatenated/ITS/cross_sample_ITS.fasta"
        )

    return outputs