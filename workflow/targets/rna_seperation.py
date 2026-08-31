from typing import Dict, List
from pathlib import Path
from snakemake.io import expand


# --------------------------
# SortMeRNA stage logic
# --------------------------

def sortmerna_stages(config: dict) -> List[str]:
    base = config["RNA"]["sortmerna_staged"]["pipeline"]

    if not base:
        raise ValueError(
            "RNA.sortmerna_staged.pipeline cannot be empty"
        )

    return [
        "_".join(base[:i + 1])
        for i in range(len(base))
    ]


def previous_stage_map(
    stages: List[str],
) -> Dict[str, str]:
    return {
        stages[0]: "decontamination",
        **{
            stages[i]: stages[i - 1]
            for i in range(1, len(stages))
        },
    }


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
    return config["databases"][f"sortmeRNA_{database_name}"]


def get_sortmerna_db_idx(
    stage: str,
    config: dict,
) -> str:
    database_name = stage.split("_")[-1]
    return config["databases"][f"sortmeRNA_{database_name}_idx"]


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
            f"{results_dir}/classified/{{sample}}/"
            f"rRNA/{{sample}}_rRNA_{read}.fastq.gz",
            sample=samples,
        )

        outputs += expand(
            f"{results_dir}/classified/{{sample}}/"
            f"non_rRNA/{{sample}}_non_rRNA_{read}.fastq.gz",
            sample=samples,
        )

    return outputs


# --------------------------
# SortMeRNA classified outputs
# --------------------------

def sortmerna_outputs(
    results_dir: Path,
    samples: List[str],
    refinement: str,
) -> List[str]:

    if refinement != "staged":
        raise ValueError(
            "Only staged SortMeRNA is currently supported "
            "by sortmerna_outputs()."
        )

    outputs: List[str] = []

    for classification in [
        "SSU",
        "LSU",
        "non_SSU",
    ]:
        for read in [1, 2]:
            outputs += expand(
                f"{results_dir}/classified/{{sample}}/"
                f"{classification}/"
                f"{{sample}}_{classification}_{read}.fastq.gz",
                sample=samples,
            )

    return outputs


# --------------------------
# RiboDetector outputs
# --------------------------

def ribodetector_outputs(
    results_dir: Path,
    samples: List[str],
    refinement: str,
) -> List[str]:

    outputs: List[str] = []

    for read in [1, 2]:
        outputs += expand(
            f"{results_dir}/intermediate/{{sample}}/"
            f"ribodetector/{{sample}}_rRNA_{read}.fastq.gz",
            sample=samples,
        )

        outputs += expand(
            f"{results_dir}/intermediate/{{sample}}/"
            f"ribodetector/{{sample}}_nonrRNA_{read}.fastq.gz",
            sample=samples,
        )

    if refinement == "bbduk":
        for classification in [
            "SSU",
            "non_SSU",
        ]:
            for read in [1, 2]:
                outputs += expand(
                    f"{results_dir}/classified/{{sample}}/"
                    f"{classification}/"
                    f"{{sample}}_{classification}_{read}.fastq.gz",
                    sample=samples,
                )

    return outputs


# --------------------------
# BBMap outputs
# --------------------------

def bbmap_outputs(
    results_dir: Path,
    samples: List[str],
) -> List[str]:

    return expand(
        f"{results_dir}/intermediate/{{sample}}/"
        f"bbmap/{{sample}}_all_reads.bam",
        sample=samples,
    )


# --------------------------
# RNA-separation targets
# --------------------------

def rna_outputs(
    results_dir: Path,
    samples: List[str],
    config: dict,
) -> List[str]:

    method = config["RNA"]["method"]
    refinement = config["RNA"]["refinement"]

    print(
        f"[RNA separation] "
        f"method={method}, refinement={refinement}"
    )

    outputs: List[str] = []

    if method == "sortmerna":
        outputs += sortmerna_outputs(
            results_dir,
            samples,
            refinement,
        )

    elif method == "ribodetector":
        outputs += ribodetector_outputs(
            results_dir,
            samples,
            refinement,
        )

    elif method == "bbmap":
        outputs += bbmap_outputs(
            results_dir,
            samples,
        )

    else:
        raise ValueError(
            f"Unsupported RNA method: {method}"
        )

    outputs += filtered_outputs(
        results_dir,
        samples,
    )

    return outputs