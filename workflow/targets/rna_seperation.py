from typing import List,Dict
from pathlib import Path
from snakemake.io import expand

# --------------------------
# sortmeRNA stage logic
# --------------------------

def sortmerna_stages(config: dict) -> List[str]:
    base = config["RNA"]["sortmerna_staged"]["pipeline"]
    return ["_".join(base[: i + 1]) for i in range(len(base))]


def previous_stage_map(stages: List[str]) -> Dict[str, str]:
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
    qc_base = results_dir / "qc" / sample
    rna_base = results_dir / "RNA" / sample

    prev = previous_stage[stage]

    if prev == "decontamination":
        return str(
            qc_base / "decontamination" / f"{sample}_R1.cleaned.fastq.gz"
        )

    return str(
        rna_base / "sortmerna" / prev / f"{sample}_{prev}.nonaligned_fwd.fq.gz"
    )


def sortmerna_input_r2(
    sample: str,
    stage: str,
    results_dir: Path,
    previous_stage: Dict[str, str],
) -> str:
    qc_base = results_dir / "qc" / sample
    rna_base = results_dir / "RNA" / sample

    prev = previous_stage[stage]

    if prev == "decontamination":
        return str(
            qc_base / "decontamination" / f"{sample}_R2.cleaned.fastq.gz"
        )

    return str(
        rna_base
        / "sortmerna"
        / prev
        / f"{sample}_{prev}.nonaligned_rev.fq.gz"
    )

# --------------------------
# Database lookup
# --------------------------

def get_sortmerna_db(stage: str, config: dict) -> str:
    key = stage.split("_")[-1]
    return config["databases"][f"sortmeRNA_{key}"]


def get_sortmerna_db_idx(stage: str, config: dict) -> str:
    key = stage.split("_")[-1]
    return config["databases"][f"sortmeRNA_{key}_idx"]

# --------------------------
# Shared outputs for the filtered rRNA and non-rRNA reads, which are common across all methods
# --------------------------
def filtered_outputs(results_dir: Path, samples: List[str]) -> List[str]:
    """
    Outputs shared across RNA separation methods.
    """

    outputs: List[str] = []

    outputs += expand(
    f"{results_dir}/classified/{{sample}}/rRNA/{{sample}}_rRNA_1.fastq.gz",
    sample=samples,
    )

    outputs += expand(
        f"{results_dir}/classified/{{sample}}/non_rRNA/{{sample}}_non_rRNA_1.fastq.gz",
        sample=samples,
    )

    return outputs

# --------------------------
# SortMeRNA outputs
# --------------------------
def sortmerna_outputs(
    results_dir: Path,
    samples: List[str],
    refinement: str,
    final_stages: List[str],
) -> List[str]:

    outputs: List[str] = []

    if refinement == "staged":
        outputs += expand(
            f"{results_dir}/RNA/{{sample}}/sortmerna/{{stage}}/{{sample}}_{{stage}}.aligned_fwd.fq.gz",
            zip,
            sample=samples,
            stage=final_stages,
        )

    elif refinement == "combined":
        outputs += expand(
            f"{results_dir}/RNA/{{sample}}/sortmerna/aligned/{{sample}}_ssu_lsu.aligned.sam",
            sample=samples,
        )

        outputs += expand(
            f"{results_dir}/RNA/{{sample}}/sortmerna/aligned/{{sample}}_SSU_filtered_taxonomy.tsv",
            sample=samples,
        )

    else:
        raise ValueError(f"Unsupported sortmerna refinement: {refinement}")

    return outputs


# --------------------------
# Ribodetector outputs
# --------------------------
def ribodetector_outputs(
    results_dir: Path,
    samples: List[str],
    refinement: str,
) -> List[str]:

    outputs: List[str] = []

    # Intermediate RiboDetector outputs

    outputs += expand(
        f"{results_dir}/intermediate/{{sample}}/ribodetector/{{sample}}_rRNA_1.fastq.gz",
        sample=samples,
    )

    outputs += expand(
        f"{results_dir}/intermediate/{{sample}}/ribodetector/{{sample}}_nonrRNA_1.fastq.gz",
        sample=samples,
    )

    # BBDuk refinement

    if refinement == "bbduk":

        outputs += expand(
            f"{results_dir}/classified/{{sample}}/SSU/{{sample}}_SSU_1.fastq.gz",
            sample=samples,
        )

        outputs += expand(
            f"{results_dir}/classified/{{sample}}/non_SSU/{{sample}}_non_SSU_1.fastq.gz",
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

    outputs: List[str] = []

    outputs += expand(
        f"{results_dir}/intermediate/{{sample}}/bbmap/{{sample}}_all_reads.bam",
        sample=samples,
    )

    return outputs

# --------------------------
# output used in Snakefile
# --------------------------
def rna_outputs(
    results_dir: Path,
    samples: List[str],
    config: dict,
    final_stages: List[str],
) -> List[str]:
    """
    Return RNA separation outputs based on configuration.
    """

    method: str = config["RNA"]["method"]
    refinement: str = config["RNA"]["refinement"]

    print(f"[RNA separation] method={method}, refinement={refinement}")

    outputs: List[str] = []

    if method == "sortmerna":
        outputs += sortmerna_outputs(
            results_dir, 
            samples, 
            refinement, 
            final_stages
        )

    elif method == "ribodetector":
        outputs += ribodetector_outputs(
            results_dir,
            samples,
            refinement,
        )

    elif method == "bbmap":
        outputs += bbmap_outputs(results_dir, samples)

    else:
        raise ValueError(f"Unsupported RNA method: {method}")

    # those refinement methods which simply add links to the original files, e.g. bbduk, will not have 
    # additional outputs, but the filtered outputs will be the same as for the other methods, 
    # so we can add those as shared across all the methods
    outputs += filtered_outputs(
        results_dir, 
        samples
        )

    return outputs