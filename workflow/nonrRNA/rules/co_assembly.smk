rule rnaSpades:
    conda:
        "../envs/spades.yaml"
    message:
        "[RNA spades] assembly across all samples"
    input:
        concatenated_fastq_r1 = f"{RESULTS_DIR}/nonrRNA/concatenated/nonrRNA_1.fastq.gz"
        concatenated_fastq_r2 = f"{RESULTS_DIR}/nonrRNA/concatenated/nonrRNA_2.fastq.gz"
    output:
        assembled_fq=f"{RESULTS_DIR}/nonrRNA/assembly/nonrRNA_assembly.fasta"
    log:
        f"{RESULTS_DIR}/nonrRNA/assembly/logs/nonrRNA_assembly.log"
    benchmark:
        f"{RESULTS_DIR}/nonrRNA/assembly/benchmarks/nonrRNA_assembly.txt"
    threads:
        config["nonrRNA"]["threads"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.genefamilies})

        spades.py -1 {input.concatenated_fastq_r1} \
            -2 {input.concatenated_fastq_r2} \
            --threads {threads} \
            --rna -o $outdir > {log} 2>&1
        """