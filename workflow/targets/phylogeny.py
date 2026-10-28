from typing import List
from pathlib import Path
from snakemake.io import expand


def phylogeny_outputs(
    results_dir: Path,
    active_domains: List[str],
    domains_with_hits: List[str],
) -> List[str]:
    """Return one main output per phylogeny rule."""

    outputs: List[str] = []

    # Domain availability.
    outputs.append(f"{results_dir}/Phylogeny/domain_status.tsv")

    # BLAST search.
    outputs += expand(
        f"{results_dir}/Phylogeny/{{domain}}/blastn/SSU_blastn_{{domain}}.tsv",
        domain=active_domains,
    )

    # Reference selection.
    outputs += expand(
        f"{results_dir}/Phylogeny/{{domain}}/blastn_postprocessing/SSU_blastn_{{domain}}_reference_candidates.tsv",
        domain=active_domains,
    )

    # Filtered BLAST availability.
    outputs += expand(
        f"{results_dir}/Phylogeny/{{domain}}/blastn_postprocessing/SSU_blastn_{{domain}}_status.txt",
        domain=active_domains,
    )

    # Reference accessions.
    outputs += expand(
        f"{results_dir}/Phylogeny/{{domain}}/references/samples_SSU_{{domain}}_accessions.txt",
        domain=domains_with_hits,
    )

    # Reference extraction.
    outputs += expand(
        f"{results_dir}/Phylogeny/{{domain}}/references/samples_SSU_{{domain}}_references.fasta",
        domain=domains_with_hits,
    )

    # Header renaming.
    outputs += expand(
        f"{results_dir}/Phylogeny/{{domain}}/references/samples_SSU_{{domain}}_renamed.fasta",
        domain=domains_with_hits,
    )

    # Combined queries and references.
    outputs += expand(
        f"{results_dir}/Phylogeny/{{domain}}/references/samples_SSU_{{domain}}_with_references.fasta",
        domain=domains_with_hits,
    )

    # Alignment.
    outputs += expand(
        f"{results_dir}/Phylogeny/{{domain}}/mafft/samples_SSU_{{domain}}_references_mafft.fasta",
        domain=domains_with_hits,
    )

    # Trimmed alignment.
    outputs += expand(
        f"{results_dir}/Phylogeny/{{domain}}/mafft/samples_SSU_{{domain}}_references_mafft_trim.fasta",
        domain=domains_with_hits,
    )

    # Tree.
    outputs += expand(
        f"{results_dir}/Phylogeny/{{domain}}/tree/samples_SSU_{{domain}}.treefile",
        domain=domains_with_hits,
    )

    return outputs