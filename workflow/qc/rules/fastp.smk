rule fastp:
    conda:
        "../envs/fastp.yaml"
    message:
        "[FastP] nucleotide quality trimming {wildcards.sample}"
    input:
        r1=raw_r1,
        r2=raw_r2,
        outdir=outdir_for_sample
    output:
        r1="{outdir}/{sample}/QC/trimmed/{sample}_R1.fastq.gz",
        r2="{outdir}/{sample}/QC/trimmed/{sample}_R2.fastq.gz",
        json="{outdir}/{sample}/QC/trimmed/{sample}_fastp.json",
        html="{outdir}/{sample}/QC/trimmed/{sample}_fastp.html",
    log:
        stdout="{outdir}/{sample}/logs/fastp.log"
    benchmark:
        "{outdir}/{sample}/benchmarks/fastp.txt"
    wildcard_constraints:
        outdir=".+"   # still correct
    params:
        options=config["qc"]["fastp"]["options"]
    threads:
        config["qc"]["fastp"].get("threads", 2)
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.r1})

        fastp \
          -i {input.r1} \
          -I {input.r2} \
          -o {output.r1} \
          -O {output.r2} \
          --json {output.json} \
          --html {output.html} \
          --thread {threads} \
          {params.options} > {log.stdout} 2>&1

        """