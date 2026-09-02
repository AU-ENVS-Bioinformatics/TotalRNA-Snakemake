from typing import Dict, List
from pathlib import Path
from snakemake.io import expand


# --------------------------
# SortMeRNA stage logic
# --------------------------

def sortmerna_stages(
    config: dict,
) -> List[str]:

    pipeline = config["RNA"]["sortmerna_staged"]["pipeline"]

    if not pipeline:
        raise ValueError(
            "RNA.sortmerna_staged.pipeline cannot be empty"
        )

    return [
        "_".join(pipeline[:i + 1])
        for i in range(len(pipeline))
    ]

def previous_stage_map(
    stages: List[str],
) -> Dict[str, str]:

    if not stages:
        raise ValueError(
            "SortMeRNA stages cannot be empty"
        )

    return {
        stages[0]: "decontamination",
        **{
            stages[i]: stages[i - 1]
            for i in range(1, len(stages))
        },
    }

# =========================================================
# SortMeRNA staged inputs
# =========================================================

def sortmerna_input_r1(
    sample: str,
    stage: str,
    results_dir: Path,
    previous_stage: Dict[str, str],
) -> str:

    previous = previous_stage[stage]

    if previous == "decontamination":
        return str(
            results_dir
            / "Quality_Control"
            / sample
            / "decontamination"
            / f"{sample}_R1.cleaned.fastq.gz"
        )

    return str(
        results_dir
        / "RNA_Classification"
        / "intermediate"
        / sample
        / "sortmerna"
        / previous
        / f"{sample}_{previous}.nonaligned_fwd.fq.gz"
    )


def sortmerna_input_r2(
    sample: str,
    stage: str,
    results_dir: Path,
    previous_stage: Dict[str, str],
) -> str:

    previous = previous_stage[stage]

    if previous == "decontamination":
        return str(
            results_dir
            / "Quality_Control"
            / sample
            / "decontamination"
            / f"{sample}_R2.cleaned.fastq.gz"
        )

    return str(
        results_dir
        / "RNA_Classification"
        / "intermediate"
        / sample
        / "sortmerna"
        / previous
        / f"{sample}_{previous}.nonaligned_rev.fq.gz"
    )


# --------------------------
# SortMeRNA database lookup
# --------------------------

def get_sortmerna_db(
    stage: str,
    config: dict,
) -> str:

    database_name = stage.split("_")[-1]
    database_key = f"sortmeRNA_{database_name}"

    if database_key not in config["databases"]:
        raise ValueError(
            f"Missing SortMeRNA database in config: {database_key}"
        )

    return config["databases"][database_key]

def get_sortmerna_db_idx(
    stage: str,
    config: dict,
) -> str:

    database_name = stage.split("_")[-1]
    database_key = f"sortmeRNA_{database_name}_idx"

    if database_key not in config["databases"]:
        raise ValueError(
            f"Missing SortMeRNA index in config: {database_key}"
        )

    return config["databases"][database_key]

# --------------------------
# Shared classified outputs
# --------------------------

def filtered_outputs(
    results_dir: Path,
    samples: List[str],
) -> List[str]:

    outputs: List[str] = []

    for read in [1, 2]:
        outputs += expand(
            f"{results_dir}/classified/{{sample}}/SSU/{{sample}}_SSU_{read}.fastq.gz",
            sample=samples,
        )

        outputs += expand(
            f"{results_dir}/classified/{{sample}}/non_rRNA/{{sample}}_non_rRNA_{read}.fastq.gz",
            sample=samples,
        )

    return outputs

# --------------------------
# ITS classified outputs
# --------------------------

def ITS_candidates(
    results_dir: Path,
    samples: List[str],
) -> List[str]:

    outputs: List[str] = []

    for read in [1, 2]:
        outputs += expand(
            f"{results_dir}/intermediate/{{sample}}/ITS_candidates/{{sample}}_ITS_candidates_{read}.fastq.gz",
            sample=samples,
        )

    return outputs

# =========================================================
# Validate RNA-separation configuration
# =========================================================

def validate_rna_config(
    config: dict,
) -> None:

    method = config["RNA"]["method"]
    refinement = config["RNA"]["refinement"]

    if method == "sortmerna":
        supported_refinements = [
            "staged",
            "combined",
        ]

        if refinement not in supported_refinements:
            raise ValueError(
                "Unsupported SortMeRNA refinement: "
                f"{refinement}. Supported refinements are: "
                + ", ".join(supported_refinements)
            )

    elif method == "ribodetector":
        if refinement != "bbduk":
            raise ValueError(
                "RiboDetector must use refinement: bbduk "
                "to produce the standardized SSU output."
            )

    elif method == "bbmap":
        # BBMap must have downstream rules that publish
        # both standardized SSU and non-rRNA FASTQs.
        pass

    elif method == "ITS":
        # ITS candidates must have downstream rules that publish
        # merged ITS candidate FASTAs.
        pass

    else:
        raise ValueError(
            f"Unsupported RNA-separation method: {method}"
        )

# =========================================================
# RNA-separation outputs used by the Snakefile
# =========================================================

def rna_outputs(
    results_dir: Path,
    samples: List[str],
    config: dict,
) -> List[str]:

    method = config["RNA"]["method"]
    refinement = config["RNA"]["refinement"]
    secondary_method = config["RNA"]["secondary_method"]
    print(f"[RNA separation] method={method}, refinement={refinement}")

    validate_rna_config(config)

    if method not in ["sortmerna", "ribodetector", "bbmap"]:
        raise ValueError(
            f"Unsupported RNA-separation method: {method}"
        )

    final_output = filtered_outputs(results_dir,samples)

    if secondary_method in ["ITS"]:
        final_output += ITS_candidates(results_dir, samples)

    return final_output