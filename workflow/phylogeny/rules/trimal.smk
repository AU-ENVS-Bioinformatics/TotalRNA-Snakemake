################################################################################
# TRIM THE SSU ALIGNMENT
################################################################################

rule trimal_ssu:
    conda:
        "../envs/trimal.yaml"
    wildcard_constraints:
        domain=PHYLOGENY_DOMAIN_PATTERN
    message:
        "[trimAl] Trimming the SSU alignment for {wildcards.domain}"
    input:
        alignment=f"{PHYLOGENY_DIR}/{{domain}}/mafft/samples_SSU_{{domain}}_references_mafft.fasta"
    output:
        trimmed=f"{PHYLOGENY_DIR}/{{domain}}/mafft/samples_SSU_{{domain}}_references_mafft_trim.fasta"
    log:
        stdout=f"{PHYLOGENY_DIR}/{{domain}}/mafft/logs/trimal.log"
    benchmark:
        f"{PHYLOGENY_DIR}/{{domain}}/mafft/benchmarks/trimal.txt"
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.trimmed}) $(dirname {log.stdout})

        trimal \
            -in {input.alignment} \
            -out {output.trimmed} \
            -automated1 \
            > {log.stdout} 2>&1
        """