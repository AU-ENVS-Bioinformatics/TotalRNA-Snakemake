rule quast:
    input:
        fasta="results/metarib/final_contigs.fasta",
        bam=expand("results/rRNA/bwa/{sample}_sorted.bam", sample=unique_samples),
    output:
        report_txt="results/qc/quast/report.txt",
        report_tsv="results/qc/quast/report.tsv",
        report_html="results/qc/quast/report.html",
        outdir=directory("results/qc/quast"),
    log:
        "logs/quast/quast.log",
    benchmark:
        "results/benchmarks/quast.txt",
    conda:
        "../envs/quast.yaml"
    params:
        extra=" ".join(config.get("quast", "")),
        comma_bam=lambda wildcards, input: ",".join(input.bam),
        space_fasta=lambda wildcards, input: " ".join([input.fasta] * len(input.bam)),
    threads: config["threads"]["quast"]
    shell:
        "quast {params.extra} "
        "--threads {threads} "
        "--bam {params.comma_bam} "
        "-o {output.outdir} "
        "{params.space_fasta} "
        "> {log} 2>&1 "
