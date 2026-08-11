################################################################################
# 9. Maximum-likelihood phylogenetic tree with IQ-TREE2
################################################################################

rule iqtree_ssu:
    conda:
        "../envs/iqtree.yaml"
    message:
        "[IQ-TREE2] construct SSU phylogenetic tree for {wildcards.sample}"
    input:
        alignment=f"{PHYLOGENY_DIR}/{{sample}}/Alignment/{{sample}}_SSU_mafft_trimmed.fasta"
    output:
        tree=f"{PHYLOGENY_DIR}/{{sample}}/IQTREE/{{sample}}_SSU.treefile",
        iqtree=f"{PHYLOGENY_DIR}/{{sample}}/IQTREE/{{sample}}_SSU.iqtree",
        model=f"{PHYLOGENY_DIR}/{{sample}}/IQTREE/{{sample}}_SSU.model.gz"
    log:
        f"{PHYLOGENY_DIR}/{{sample}}/logs/{{sample}}_iqtree.log"
    benchmark:
        f"{PHYLOGENY_DIR}/{{sample}}/benchmarks/{{sample}}_iqtree.txt"
    params:
        prefix=f"{PHYLOGENY_DIR}/{{sample}}/IQTREE/{{sample}}_SSU",
        options=config["phylogeny"]["iqtree"]["options"]
    threads:
        config["phylogeny"]["iqtree"].get("threads", 8)
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.tree})

        iqtree2 \
            -s {input.alignment} \
            -pre {params.prefix} \
            -T {threads} \
            > {log} 2>&1
        """