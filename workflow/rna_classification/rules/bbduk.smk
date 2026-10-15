rule bbduk_ssu:
    conda:
        "../envs/bbmap.yaml"
    message:
        "[BBDuk] Extracting SSU reads for {wildcards.sample}"
    input:
        rRNA_r1=rules.ribodetector.output.rna_r1,
        rRNA_r2=rules.ribodetector.output.rna_r2,
        ref=config["databases"]["sortmeRNA_ssu"]
    output:
        rRNA_ssu_r1=f"{RNA_CLASSIFIED_DIR}/{{sample}}/SSU/{{sample}}_SSU_1.fastq.gz",
        rRNA_ssu_r2=f"{RNA_CLASSIFIED_DIR}/{{sample}}/SSU/{{sample}}_SSU_2.fastq.gz",
        rRNA_non_ssu_r1=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/ribodetector/non_SSU/{{sample}}_non_SSU_1.fastq.gz",
        rRNA_non_ssu_r2=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/ribodetector/non_SSU/{{sample}}_non_SSU_2.fastq.gz",
    log:
        stdout=f"{RNA_CLASSIFIED_DIR}/{{sample}}/logs/bbduk_ssu.log"
    benchmark:
        f"{RNA_CLASSIFIED_DIR}/{{sample}}/benchmarks/bbduk_ssu.txt"
    params:
        options=config["RNA"]["bbduk"]["options"]
    threads:
        config["RNA"]["bbduk"]["threads"]
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.rRNA_ssu_r1})

        bbduk.sh \
            in={input.rRNA_r1} \
            in2={input.rRNA_r2} \
            outm={output.rRNA_ssu_r1} \
            outm2={output.rRNA_ssu_r2} \
            out={output.rRNA_non_ssu_r1} \
            out2={output.rRNA_non_ssu_r2} \
            ref={input.ref} \
            threads={threads} \
            {params.options} \
            > {log.stdout} 2>&1
        """
