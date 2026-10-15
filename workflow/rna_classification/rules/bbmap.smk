rule bbmap:
    conda:
        "../envs/bbmap.yaml"
    message:
        "[BBMap] Running BBMap for {wildcards.sample}"
    input:
        cleaned_r1=f"{QC_DIR}/{{sample}}/decontamination/{{sample}}_R1.cleaned.fastq.gz",
        cleaned_r2=f"{QC_DIR}/{{sample}}/decontamination/{{sample}}_R2.cleaned.fastq.gz",
        ref=config["databases"]["sortmeRNA_ssu_lsu"]
    output:
        all_RNA_reads_bam=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/bbmap/RNA/{{sample}}_all_reads.bam",
    log:
        stdout=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/logs/bbmap_RNA.log"
    benchmark:
        f"{RNA_INTERMEDIATE_DIR}/{{sample}}/benchmarks/bbmap_RNA.txt"
    params:
        options=config["RNA"]["bbmap"]["options"]
    threads:
        config["RNA"]["bbmap"]["threads"]
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.all_RNA_reads_bam})
        mkdir -p $(dirname {log.stdout})

        bbmap.sh \
            in={input.cleaned_r1} \
            in2={input.cleaned_r2} \
            ref={input.ref} \
            out={output.all_RNA_reads_bam} \
            threads={threads} \
            {params.options} \
            > {log.stdout} 2>&1
        """
