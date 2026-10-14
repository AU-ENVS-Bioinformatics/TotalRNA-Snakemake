################################################################################
# 9. Maximum-likelihood phylogenetic tree with IQ-TREE2
################################################################################
rule iqtree_ssu:
    conda:
        "../envs/iqtree.yaml"
    wildcard_constraints:
        kingdom=BARRNAP_KINGDOM_PATTERN
    message:
        "[IQ-TREE2] Construct the cross-sample SSU phylogenetic tree"
    input:
        trimmed=f"{PHYLOGENY_DIR}/{{kingdom}}/mafft/samples_SSU_{{kingdom}}_references_mafft_trim.fasta"
    output:
        tree=f"{PHYLOGENY_DIR}/{{kingdom}}/tree/samples_SSU_{{kingdom}}.treefile",
        iqtree=f"{PHYLOGENY_DIR}/{{kingdom}}/tree/samples_SSU_{{kingdom}}.iqtree",
        model=f"{PHYLOGENY_DIR}/{{kingdom}}/tree/samples_SSU_{{kingdom}}.model.gz"
    log:
        stdout=f"{PHYLOGENY_DIR}/{{kingdom}}/tree/logs/iqtree.log",
        time=f"{PHYLOGENY_DIR}/{{kingdom}}/tree/benchmarks/iqtree_time.txt"
    benchmark:
        f"{PHYLOGENY_DIR}/{{kingdom}}/tree/benchmarks/iqtree.txt"
    params:
        prefix=f"{PHYLOGENY_DIR}/{{kingdom}}/tree/samples_SSU_{{kingdom}}",
        options=lambda wc: config["phylogeny"]["iqtree"][wc.kingdom].get("options"),
        seed=config["phylogeny"]["iqtree"].get("seed", 42)
    threads:
        config["phylogeny"]["iqtree"].get("threads", 8)
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.tree})

        /usr/bin/time -v -o {log.time} \
            iqtree2 \
                -s {input.trimmed} \
                -pre {params.prefix} \
                -T {threads} \
                {params.options} \
                -seed {params.seed} \
                > {log.stdout} 2>&1
        """