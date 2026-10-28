################################################################################
# MAXIMUM-LIKELIHOOD SSU TREE
################################################################################

rule iqtree_ssu:
    conda:
        "../envs/iqtree.yaml"
    wildcard_constraints:
        domain=PHYLOGENY_DOMAIN_PATTERN
    message:
        "[IQ-TREE2] Constructing the SSU tree for {wildcards.domain}"
    input:
        trimmed=f"{PHYLOGENY_DIR}/{{domain}}/mafft/samples_SSU_{{domain}}_references_mafft_trim.fasta"
    output:
        tree=f"{PHYLOGENY_DIR}/{{domain}}/tree/samples_SSU_{{domain}}.treefile",
        iqtree=f"{PHYLOGENY_DIR}/{{domain}}/tree/samples_SSU_{{domain}}.iqtree"
    log:
        stdout=f"{PHYLOGENY_DIR}/{{domain}}/tree/logs/iqtree.log",
        time=f"{PHYLOGENY_DIR}/{{domain}}/tree/benchmarks/iqtree_time.txt"
    benchmark:
        f"{PHYLOGENY_DIR}/{{domain}}/tree/benchmarks/iqtree.txt"
    params:
        prefix=f"{PHYLOGENY_DIR}/{{domain}}/tree/samples_SSU_{{domain}}",
        options=lambda wc: config["phylogeny"]["iqtree"][wc.domain].get("options", ""),
        seed=config["phylogeny"]["iqtree"].get("seed", 42)
    threads:
        config["phylogeny"]["iqtree"].get("threads", 12)
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.tree}) $(dirname {log.stdout}) $(dirname {log.time})

        /usr/bin/time -v -o {log.time} \
            iqtree2 \
                -s {input.trimmed} \
                -pre {params.prefix} \
                -T {threads} \
                {params.options} \
                -seed {params.seed} \
                > {log.stdout} 2>&1
        """