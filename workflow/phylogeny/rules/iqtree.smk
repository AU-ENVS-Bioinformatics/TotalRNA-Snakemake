################################################################################
# 9. Maximum-likelihood phylogenetic tree with IQ-TREE2
################################################################################

rule iqtree_ssu:
    conda:
        "../envs/iqtree.yaml"
    message:
        "[IQ-TREE2] construct SSU phylogenetic tree across samples and SSU extracted references"
    input:
        trimmed=f"{PHYLOGENY_DIR}/mafft/cross_sample_SSU_reference_mafft_trim.fasta"
    output:
        tree=f"{PHYLOGENY_DIR}/tree/cross_sample_SSU.treefile",
        iqtree=f"{PHYLOGENY_DIR}/tree/cross_sample_SSU.iqtree",
        model=f"{PHYLOGENY_DIR}/tree/cross_sample_SSU.model.gz"
    log:
        stdout=f"{PHYLOGENY_DIR}/tree/logs/iqtree.log",
        benchmark_file=f"{PHYLOGENY_DIR}/tree/benchmarks/iqtree_time.txt",
    benchmark:
        f"{PHYLOGENY_DIR}/tree/benchmarks/iqtree.txt"
    params:
        prefix=f"{PHYLOGENY_DIR}/tree/cross_sample_SSU",
        options=config["phylogeny"]["iqtree"]["options"]
    threads:
        config["phylogeny"]["iqtree"].get("threads", 8)
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.tree})

        /usr/bin/time -v -o {log.benchmark_file} \
            iqtree2 \
                -s {input.trimmed} \
                -pre {params.prefix} \
                -T {threads} \
                > {log.stdout} 2>&1
        """