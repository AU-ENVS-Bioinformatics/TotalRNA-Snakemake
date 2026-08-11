################################################################################
# 7. Multiple sequence alignment with MAFFT
################################################################################

rule mafft_ssu:
    conda:
        "../envs/mafft.yaml"
    message:
        "[MAFFT] align reconstructed SSUs and SILVA references for {wildcards.sample}"
    input:
        fasta=f"{PHYLOGENY_DIR}/{{sample}}/MAS/{{sample}}_SSU_with_references.fasta"
    output:
        alignment=f"{PHYLOGENY_DIR}/{{sample}}/MAS/{{sample}}_SSU_mafft.fasta"
    log:
        f"{PHYLOGENY_DIR}/{{sample}}/logs/{{sample}}_mafft.log"
    benchmark:
        f"{PHYLOGENY_DIR}/{{sample}}/benchmarks/{{sample}}_mafft.txt"
    threads:
        config["phylogeny"]["mafft"].get("threads", 8)
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.alignment})

        mafft \
            --auto \
            --thread {threads} \
            {input.fasta} \
            > {output.alignment} \
            2> {log}
        """

################################################################################
# 8. Trim alignment with trimAl
################################################################################

rule trimal_ssu:
    conda:
        "../envs/trimal.yaml"
    message:
        "[trimAl] trim SSU alignment for {wildcards.sample}"
    input:
        alignment=f"{PHYLOGENY_DIR}/{{sample}}/Alignment/{{sample}}_SSU_mafft.fasta"
    output:
        trimmed=f"{PHYLOGENY_DIR}/{{sample}}/Alignment/{{sample}}_SSU_mafft_trimmed.fasta"
    log:
        f"{PHYLOGENY_DIR}/{{sample}}/logs/{{sample}}_trimal.log"
    benchmark:
        f"{PHYLOGENY_DIR}/{{sample}}/benchmarks/{{sample}}_trimal.txt"
    shell:
        r"""
        set -euo pipefail

        trimal \
            -in {input.alignment} \
            -out {output.trimmed} \
            -automated1 \
            > {log} 2>&1
        """
