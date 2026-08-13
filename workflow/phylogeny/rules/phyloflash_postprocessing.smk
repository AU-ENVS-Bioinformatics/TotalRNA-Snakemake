################################################################################
# 2. PREFIX PHYLOFLASH SEQUENCE IDS WITH SAMPLE NAME
################################################################################

rule prefix_phyloflash_ids:
    message:
        "[PhyloFlash postprocessing] prefix reconstructed SSU IDs for {wildcards.sample}"
    input:
        fasta=f"{PHYLOGENY_DIR}/{{sample}}/phyloflash/{{sample}}.all.final.fasta"
    output:
        fasta=f"{PHYLOGENY_DIR}/{{sample}}/cross_sample/{{sample}}_SSU_prefixed.fasta"
    params:
        sample=lambda wildcards: wildcards.sample
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname "{output.fasta}")"

        awk -v sample="{params.sample}" '
            /^>/ {{
                sub(/^>/, ">" sample "__")
            }}
            {{ print }}
        ' "{input.fasta}" > "{output.fasta}"
        """

################################################################################
# 3. COMBINE PREFIXED SSU RECONSTRUCTIONS ACROSS SAMPLES
################################################################################

rule combine_reconstructed_ssu:
    message:
        "[Phylogeny] combine reconstructed SSUs across all samples"
    input:
        fasta=expand(
            f"{PHYLOGENY_DIR}/{{sample}}/phyloflash/{{sample}}.all.final.fasta",
            sample=SAMPLES
        )
    output:
        fasta=f"{PHYLOGENY_DIR}/cross_sample/reconstructed_SSU_all_samples.fasta"
    benchmark:
        f"{PHYLOGENY_DIR}/cross_sample/benchmarks/combine_reconstructed_ssu.txt"
    shell:
        r"""
        set -euo pipefail

        mkdir -p \
            "$(dirname "{output.fasta}")"

        cat {input.fasta:q} > "{output.fasta}"
        """