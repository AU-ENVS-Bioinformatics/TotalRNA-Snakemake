################################################################################
# ALIGN SAMPLE-DERIVED SSUS AND REFERENCES
################################################################################

rule mafft_ssu:
    conda:
        "../envs/mafft.yaml"
    wildcard_constraints:
        domain=PHYLOGENY_DOMAIN_PATTERN
    message:
        "[MAFFT] Aligning SSUs and references for {wildcards.domain}"
    input:
        fasta=f"{PHYLOGENY_DIR}/{{domain}}/references/samples_SSU_{{domain}}_with_references.fasta"
    output:
        alignment=f"{PHYLOGENY_DIR}/{{domain}}/mafft/samples_SSU_{{domain}}_references_mafft.fasta"
    log:
        stdout=f"{PHYLOGENY_DIR}/{{domain}}/mafft/logs/mafft.log",
        time=f"{PHYLOGENY_DIR}/{{domain}}/mafft/benchmarks/mafft_time.txt"
    benchmark:
        f"{PHYLOGENY_DIR}/{{domain}}/mafft/benchmarks/mafft.txt"
    threads:
        config["phylogeny"]["mafft"].get("threads", 12)
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.alignment}) $(dirname {log.stdout}) $(dirname {log.time})

        /usr/bin/time -v -o {log.time} \
            mafft \
                --auto \
                --thread {threads} \
                {input.fasta} \
                > {output.alignment} \
                2> {log.stdout}
        """