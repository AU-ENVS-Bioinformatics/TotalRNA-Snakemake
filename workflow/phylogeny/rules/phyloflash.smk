################################################################################
# 1. PHYLOFLASH SSU RECONSTRUCTIONS FOR PHYLOGENETIC ANALYSIS
################################################################################

PHYLOGENY_DIR = f"{RESULTS_DIR}/phylogeny"

rule phyloflash:
    conda:
        "../envs/phyloflash.yaml"
    message:
        "[PhyloFlash] reconstruct the SSU rRNAs and explore phylogenetic composition of {wildcards.sample}"
    input:
        rRNA_ssu_r1=f"{RESULTS_DIR}/rRNA/{{sample}}/filtered/{{sample}}_rRNA_1.fastq.gz",
        rRNA_ssu_r2=f"{RESULTS_DIR}/rRNA/{{sample}}/filtered/{{sample}}_rRNA_2.fastq.gz",
    output:
        Phyloflash_fasta=f"{PHYLOGENY_DIR}/{{sample}}/phyloflash/{{sample}}.all.final.fasta",
    log:
        stdout=f"{PHYLOGENY_DIR}/{{sample}}/logs/{{sample}}_phyloflash.log",
        benchmark_file=f"{PHYLOGENY_DIR}/{{sample}}/benchmarks/{{sample}}_phyloflash_time.txt",
    benchmark:
        f"{PHYLOGENY_DIR}/{{sample}}/benchmarks/{{sample}}_phyloflash.txt"
    params:
        dir=f"{PHYLOGENY_DIR}/{{sample}}/phyloflash/",
        db=config["databases"]["phyloflash"],
        options=config["phylogeny"]["phyloflash"]["options"]
    threads:
        config["phylogeny"]["phyloflash"].get("threads", 2)
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.dir}
                
        cd {params.dir}

        /usr/bin/time -v -o {log.benchmark_file} \
            phyloFlash.pl \
                -lib {wildcards.sample} \
                -read1 {input.rRNA_ssu_r1} \
                -read2 {input.rRNA_ssu_r2} \
                -dbhome {params.db} \
                -CPUs {threads} \
                {params.options} \
                > {log.stdout} 2>&1
        """