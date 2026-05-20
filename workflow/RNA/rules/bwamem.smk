
#bbmap.sh ref=combined.fa in=reads.fq out=mapped.sam ambig=best
rule bbmap:
    conda:
        "../envs/bbmap.yaml"
    message:
        "[BBMap] Running BBMap for {wildcards.sample}"
    input:
        cleaned_r1=f"{RESULTS_DIR}/qc/{{sample}}/decontamination/{{sample}}_R1.cleaned.fastq.gz",
        cleaned_r2=f"{RESULTS_DIR}/qc/{{sample}}/decontamination/{{sample}}_R2.cleaned.fastq.gz",
        ref=config["databases"]["bbmap_rRNA_reference"]
    output:
        aligned=f"{RESULTS_DIR}/RNA/{{sample}}/bwamem/{{sample}}_aligned.sam",
    log:
        stdout = f"{RESULTS_DIR}/RNA/{{sample}}/logs/bwamem.log"
    benchmark:
        f"{RESULTS_DIR}/RNA/{{sample}}/benchmarks/bwamem.txt"
    params:
        options=config["RNA"]["bwamem"]["options"]
    threads:
        config["RNA"]["bwamem"]["threads"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.aligned})

        bwa-mem2 mem \
            -t {threads} \
            {params.genome} \
            {input.trim_r1} \
            {input.trim_r2} \
            -o {output.aligned} \
            > {log.stdout} 2>&1
        """