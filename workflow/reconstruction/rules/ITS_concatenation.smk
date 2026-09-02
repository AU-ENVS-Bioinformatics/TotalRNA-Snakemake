################################################################################
# 2. PREFIX PHYLOFLASH SEQUENCE IDS WITH SAMPLE NAME
################################################################################

rule prefix_contigs_ids:
    message:
        "[ITS postprocessing] Prefixing reconstructed ITS IDs for {wildcards.sample}"
    input:
        fasta=f"{ASSEMBLY_DIR}/{{sample}}/ITS_candidates/ITSx/{{sample}}_ITSx.full.fasta"
    output:
        fasta=f"{ASSEMBLY_DIR}/{{sample}}/ITS_candidates/ITSx/{{sample}}_ITS_prefixed.fasta"
    params:
        sample=lambda wildcards: wildcards.sample
    shell:
        r"""
        set -euo pipefail

        awk -v sample="{params.sample}" '
            /^>/ {{
                sub(/^>/, ">" sample "__")
            }}
            {{ print }}
        ' "{input.fasta}" > "{output.fasta}"
        """

################################################################################
# 3. COMBINE PREFIXED ITS RECONSTRUCTIONS ACROSS SAMPLES
################################################################################

rule combine_reconstructed_its:
    message:
        "[ITS postprocessing] Combining reconstructed ITS across all samples"
    input:
        fasta=expand(
            f"{ASSEMBLY_DIR}/{{sample}}/ITS_candidates/ITSx/{{sample}}_ITS_prefixed.fasta",
            sample=SAMPLES
        )
    output:
        fasta=f"{ASSEMBLY_DIR}/concatenated/ITS/cross_sample_ITS.fasta"
    shell:
        r"""
        set -euo pipefail

        mkdir -p \
            "$(dirname "{output.fasta}")"

        cat {input.fasta:q} > "{output.fasta}"
        """