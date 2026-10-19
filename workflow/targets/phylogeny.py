from typing import List
from pathlib import Path
from snakemake.io import expand

def phylogeny_outputs(
    results_dir: Path,
    barrnap_kingdoms: List[str],
    kingdoms_with_blast_hits: List[str],
) -> List[str]:
    """
    Return final expected SSU phylogeny outputs.
    """

    outputs: List[str] = []

    #step 1 blastn of reconstructed sequence for the kingdoms bac, arc, fun
    outputs += expand(
        f"{results_dir}/Phylogeny/{{kingdom}}/blastn/SSU_blastn_{{kingdom}}.tsv",
        kingdom=barrnap_kingdoms,
    )

    #Step 2 blastn postprocessing of reconstructed sequence for the kingdoms bac, arc, fun
    outputs += expand(
        f"{results_dir}/Phylogeny/{{kingdom}}/blastn_postprocessing/SSU_blastn_{{kingdom}}_reference_candidates.tsv",
        kingdom=barrnap_kingdoms,
    )

    #Step 3 extract reference IDs of top blast hits for the kingdoms bac, arc, fun
    outputs += expand(
        f"{results_dir}/Phylogeny/{{kingdom}}/references/samples_SSU_{{kingdom}}_accessions.txt",
        kingdom=kingdoms_with_blast_hits,
    )

    #Step 4-5 Extract and rename reference genoms based on the accesion id's for the kingdoms bac, arc, fun
    outputs += expand(
        f"{results_dir}/Phylogeny/{{kingdom}}/references/samples_SSU_{{kingdom}}_references.fasta",
        kingdom=kingdoms_with_blast_hits,
    )

    outputs += expand(
        f"{results_dir}/Phylogeny/{{kingdom}}/references/samples_SSU_{{kingdom}}_renamed.fasta",
        kingdom=kingdoms_with_blast_hits,
    )

    #Step 6 Combine extracted top reference genomes witgh reconstructed SSU sequences
    outputs += expand(
        f"{results_dir}/Phylogeny/{{kingdom}}/references/samples_SSU_{{kingdom}}_with_references.fasta",
        kingdom=kingdoms_with_blast_hits,
    )

    #Step 7-8 Multiple sequence alignment of the combined sequences with MAFFT followed by trimming
    outputs += expand(
        f"{results_dir}/Phylogeny/{{kingdom}}/mafft/samples_SSU_{{kingdom}}_references_mafft.fasta",
        kingdom=kingdoms_with_blast_hits,
    )

    outputs += expand(
        f"{results_dir}/Phylogeny/{{kingdom}}/mafft/samples_SSU_{{kingdom}}_references_mafft_trim.fasta",
        kingdom=kingdoms_with_blast_hits,
    )

    # Step 8 Creation of the maximum likelihood phylogenetic tree with IQ-TREE2
    outputs += expand(
        f"{results_dir}/Phylogeny/{{kingdom}}/tree/samples_SSU_{{kingdom}}.treefile",
        kingdom=kingdoms_with_blast_hits,
    )

    return outputs