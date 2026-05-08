rule sortmerna:
    input:
        r1=clean_r1,
        r2=clean_r2
    output:
        r1="{out}/{sample}/QC/rrna_filtered/{sample}_R1.fastq",
        r2="{out}/{sample}/QC/rrna_filtered/{sample}_R2.fastq"
    log:
        "{out}/{sample}/logs/sortmerna.log"
    benchmark:
        "{out}/{sample}/benchmarks/sortmerna.txt"
    params:
        db1=config["databases"]["sortmeRNA_ssu"],
        db2=config["databases"]["sortmeRNA_lsu"],
        out=lambda wc: outdir(wc.sample)
    threads: config["qc"]["sortmerna"]["threads"]
    shell:
        """
        mkdir -p {params.out}/{wildcards.sample}/QC/rrna_filtered \
                 $(dirname {log})

        sortmerna \
          --ref {params.db1} \
          --ref {params.db2} \
          --reads {input.r1} \
          --reads {input.r2} \
          --paired_in \
          --fastx \
          --other {params.out}/{wildcards.sample}/QC/rrna_filtered/{wildcards.sample} \
          --threads {threads} \
          > {log} 2>&1
    """