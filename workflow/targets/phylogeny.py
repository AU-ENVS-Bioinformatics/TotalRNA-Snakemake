from typing import List
from pathlib import Path
from snakemake.io import expand

def phylogeny_outputs(
    results_dir: Path,
    samples: List[str],
    config: dict,
) -> List[str]:
    """
    Return expected phylogeny outputs.
    """

    method: str = config["phylogeny"]["reconstruction_method"]
    outputs: List[str] = []

    if method in ["phyloflash"]:

        # Per-sample kraken outputs
        outputs += expand(
            f"{results_dir}/phylogeny/{{sample}}/phyloflash/{{sample}}.all.final.fasta",
            sample=samples,
        )

        # Sample-prefixed SSU reconstructions for cross-sample analysis
        outputs += expand(
            f"{results_dir}/phylogeny/{{sample}}/cross_sample/{{sample}}_SSU_prefixed.fasta",
            sample=samples,
        )

        outputs += [
            f"{results_dir}/phylogeny/cross_sample/reconstructed_SSU_all_samples.fasta"
        ]

        outputs += [
                    f"{results_dir}/phylogeny/cross_sample/blastn/cross_sample_SSU_blastn.tsv"
                ]

        outputs += [
                f"{results_dir}/phylogeny/phyloflashreference/cross_sample_acc_SSU_with_references.fasta"
            ]

        outputs += [
                f"{results_dir}/phylogeny/mafft/cross_sample_SSU_reference_mafft.fasta"
            ]

        outputs += [
                f"{results_dir}/phylogeny/mafft/cross_sample_SSU_reference_mafft_trim.fasta"
            ]

        outputs += [
                f"{results_dir}/phylogeny/tree/cross_sample_SSU.treefile"
            ]

        return outputs
    else:
        raise ValueError(f"Unsupported phylogeny method: {method}")