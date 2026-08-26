################################################################################
# 1. RNA SPADES SSU RECONSTRUCTION FOR PHYLOGENETIC ANALYSIS
################################################################################

PHYLOGENY_DIR = f"{RESULTS_DIR}/phylogeny"

if config["phylogeny"]["reconstruction_method"] in ["spades","rnaspades"]:
    rule rnaspades:
        conda:
            "../envs/spades.yaml"
        message:
            "[rnaSpades] reconstruct the SSU rRNAs for {wildcards.sample}"
        input:
            rRNA_ssu_r1=f"{RESULTS_DIR}/rRNA/{{sample}}/filtered/{{sample}}_rRNA_1.fastq.gz",
            rRNA_ssu_r2=f"{RESULTS_DIR}/rRNA/{{sample}}/filtered/{{sample}}_rRNA_2.fastq.gz",
        output:
            Spades_fasta=f"{PHYLOGENY_DIR}/{{sample}}/reconstructed/transcripts.fasta",
        log:
            stdout=f"{PHYLOGENY_DIR}/{{sample}}/logs/{{sample}}_spades.log",
            benchmark_file=f"{PHYLOGENY_DIR}/{{sample}}/benchmarks/{{sample}}_spades_time.txt",
        benchmark:
            f"{PHYLOGENY_DIR}/{{sample}}/benchmarks/{{sample}}_rnaspades.txt"
        params:
            options=config["phylogeny"]["rnaspades"]["options"],
            outdir=f"{PHYLOGENY_DIR}/{{sample}}/reconstructed/"
        threads:
            config["phylogeny"]["rnaspades"].get("threads", 12)
        shell:
            r"""
            set -euo pipefail
            mkdir -p {params.outdir}

            /usr/bin/time -v -o {log.benchmark_file} \
                rnaspades.py \
                    -1 {input.rRNA_ssu_r1} \
                    -2 {input.rRNA_ssu_r2} \
                    -t {threads} \
                    -o {params.outdir} \
                    {params.options}
                    > {log.stdout} 2>&1
            """