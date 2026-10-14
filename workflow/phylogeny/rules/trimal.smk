################################################################################
# 7. Trim the multiple sequence alignment with trimAl
################################################################################

rule trimal_ssu:
    conda:
        "../envs/trimal.yaml"
    wildcard_constraints:
        kingdom=BARRNAP_KINGDOM_PATTERN
    message:
        "[trimAl] Trim the cross-sample SSU alignment"
    input:
        alignment=f"{PHYLOGENY_DIR}/{{kingdom}}/mafft/samples_SSU_{{kingdom}}_references_mafft.fasta"
    output:
        trimmed=f"{PHYLOGENY_DIR}/{{kingdom}}/mafft/samples_SSU_{{kingdom}}_references_mafft_trim.fasta"
    log:
        stdout=f"{PHYLOGENY_DIR}/{{kingdom}}/mafft/logs/trimal.log"
    benchmark:
        f"{PHYLOGENY_DIR}/{{kingdom}}/mafft/benchmarks/trimal.txt"
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.trimmed})

        trimal \
            -in {input.alignment} \
            -out {output.trimmed} \
            -automated1 \
            > "{log.stdout}" 2>&1
        """