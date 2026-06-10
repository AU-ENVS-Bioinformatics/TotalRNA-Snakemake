from typing import List
from pathlib import Path
from snakemake.io import expand


# --------------------------
# Shared outputs for the filtered rRNA and non-rRNA reads, which are common across all methods
# --------------------------
def filtered_outputs(results_dir: Path, samples: List[str]) -> List[str]:
    """
    Outputs shared across RNA separation methods.
    """

    outputs: List[str] = []

    outputs += expand(
        f"{results_dir}/rRNA/{{sample}}/filtered/{{sample}}_rRNA_1.fastq.gz",
        sample=samples,
    )

    outputs += expand(
        f"{results_dir}/nonrRNA/{{sample}}/filtered/{{sample}}_nonrRNA_1.fastq.gz",
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
) -> List[str]:

    outputs: List[str] = []

    # Only specify read 1, will be the same for read 2
    outputs += expand(
        f"{results_dir}/RNA/{{sample}}/ribodetector/{{sample}}_rRNA_1.fastq.gz",
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
        f"{results_dir}/RNA/{{sample}}/bbmap/{{sample}}_all_reads.bam",
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
            results_dir, samples, refinement, final_stages
        )

    elif method == "ribodetector":
        outputs += ribodetector_outputs(results_dir, samples)

    elif method == "bbmap":
        outputs += bbmap_outputs(results_dir, samples)

    else:
        raise ValueError(f"Unsupported RNA method: {method}")

    # those refinement methods which simply add links to the original files, e.g. bbduk, will not have 
    # additional outputs, but the filtered outputs will be the same as for the other methods, 
    # so we can add those as shared across all the methods
    outputs += filtered_outputs(results_dir, samples)

    return outputs