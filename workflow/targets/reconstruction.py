from typing import List
from pathlib import Path
from snakemake.io import expand

from utils.reconstruction import (
    validate_reconstruction_config,
    coassembly_targets,
    coassembly_fasta,
)


def reconstruction_outputs(
    results_dir: Path,
    samples: List[str],
    config: dict,
) -> List[str]:
    """Return one main output per reconstruction rule."""

    validate_reconstruction_config(config)

    assembly_refinement = config["reconstruction"].get("refinement", "coassembly")
    secondary_method = config["RNA_classification"].get("secondary_method", "none")

    outputs: List[str] = []

    # Co-assemblies and postprocessing.
    if assembly_refinement in {"coassembly", "both", "all"}:
        targets = coassembly_targets(config)

        # Assemblies.
        outputs += [
            coassembly_fasta(results_dir=results_dir, config=config, target=target)
            for target in targets
        ]

        # Length filtering.
        outputs += expand(
            f"{results_dir}/Assembly/coassembly/{{target}}/postprocessing/{{target}}_length_filtered.fasta",
            target=targets,
        )

        # Exact dereplication.
        outputs += expand(
            f"{results_dir}/Assembly/coassembly/{{target}}/postprocessing/{{target}}_dereplicated.fasta",
            target=targets,
        )

        # Rfam search.
        outputs += expand(
            f"{results_dir}/Assembly/coassembly/{{target}}/postprocessing/{{target}}_rfam_hits.tbl",
            target=targets,
        )

        # Non-rRNA filtering.
        if "non_rRNA" in targets:
            outputs.append(f"{results_dir}/Assembly/coassembly/non_rRNA/final/transcript.fasta")

        if "SSU" in targets:
            # SSU filtering.
            outputs.append(f"{results_dir}/Assembly/coassembly/SSU/postprocessing/SSU_rfam_filtered.fasta")

            # Orientation and chimera screening.
            outputs.append(f"{results_dir}/Assembly/coassembly/SSU/final/transcript.fasta")

    # Sample-specific assemblies and postprocessing.
    if assembly_refinement in {"samples", "both", "all"}:
        # SSU assemblies.
        outputs += expand(
            f"{results_dir}/Assembly/cross_assembly/{{sample}}/SSU/rnaspades/transcripts.fasta",
            sample=samples,
        )

        # Prefix sample identifiers.
        outputs += expand(
            f"{results_dir}/Assembly/cross_assembly/{{sample}}/SSU/rnaspades/transcripts_prefixed.fasta",
            sample=samples,
        )

        # Concatenation.
        outputs.append(f"{results_dir}/Assembly/cross_assembly/concatenated/SSU/cross_sample_SSU.fasta")

        # Length filtering.
        outputs.append(f"{results_dir}/Assembly/cross_assembly/SSU/postprocessing/SSU_length_filtered.fasta")

        # Exact dereplication.
        outputs.append(f"{results_dir}/Assembly/cross_assembly/SSU/postprocessing/SSU_dereplicated.fasta")

        # Rfam search.
        outputs.append(f"{results_dir}/Assembly/cross_assembly/SSU/postprocessing/SSU_rfam_hits.tbl")

        # SSU filtering.
        outputs.append(f"{results_dir}/Assembly/cross_assembly/SSU/postprocessing/SSU_rfam_filtered.fasta")

        # Domain splitting.
        outputs += expand(
            f"{results_dir}/Assembly/cross_assembly/SSU/domains/{{domain}}.fasta",
            domain=["archaea", "bacteria", "eukaryotes"],
        )

        # Domain-specific orientation and chimera screening.
        outputs += expand(
            f"{results_dir}/Assembly/cross_assembly/SSU/domains/{{domain}}/final/transcript.fasta",
            domain=["archaea", "bacteria", "eukaryotes"],
        )

        # Optional ITS-candidate processing.
        if secondary_method == "ITS":
            # ITS-candidate assemblies.
            outputs += expand(
                f"{results_dir}/Assembly/cross_assembly/{{sample}}/ITS_candidates/rnaspades/transcripts.fasta",
                sample=samples,
            )

            # ITSx extraction.
            outputs += expand(
                f"{results_dir}/Assembly/cross_assembly/{{sample}}/ITS_candidates/ITSx/{{sample}}_ITSx.full.fasta",
                sample=samples,
            )

            # Prefix sample identifiers.
            outputs += expand(
                f"{results_dir}/Assembly/cross_assembly/{{sample}}/ITS_candidates/ITSx/{{sample}}_ITSx.full_prefixed.fasta",
                sample=samples,
            )

            # Concatenation.
            outputs.append(f"{results_dir}/Assembly/cross_assembly/concatenated/ITS/cross_sample_ITS.fasta")

    return outputs