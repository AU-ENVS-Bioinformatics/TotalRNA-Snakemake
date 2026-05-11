rule kraken2:
    conda:
        "../envs/kraken2.yaml"
    message:
        "[Kraken2] estimate species composition for {wildcards.sample}"
    input:
        trim_r1="{outdir}/{sample}/QC/trimmed/{sample}_R1.fastq.gz",
        trim_r2="{outdir}/{sample}/QC/trimmed/{sample}_R2.fastq.gz",
    output:
        report="{outdir}/{sample}/QC/classification/{sample}.report",
        kraken="{outdir}/{sample}/QC/classification/{sample}.kraken"
    log:
        stdout = "{outdir}/{sample}/logs/kraken2.log"
    benchmark:
        "{outdir}/{sample}/benchmarks/kraken2.txt"
    params:
        db=config["databases"]["kraken_db"],
        options=config["qc"]["kraken2"]["options"],
    threads:
        config["qc"]["kraken2"].get("threads", 2)
    wildcard_constraints:
        outdir=".+"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.report})

        kraken2 \
          --db {params.db} \
          --threads {threads} \
          --report {output.report} \
          --output {output.kraken} \
          {params.options} \
          {input.trim_r1} {input.trim_r2} \
          > {log.stdout} 2>&1
        """