rule bbduk_its_candidates:
    conda:
        "../envs/bbmap.yaml"
    message:
        "[BBDuk] Extracting ITS candidate reads for {wildcards.sample}"
    input:
        cleaned_r1=f"{QC_DIR}/{{sample}}/decontamination/{{sample}}_R1.cleaned.fastq.gz",
        cleaned_r2=f"{QC_DIR}/{{sample}}/decontamination/{{sample}}_R2.cleaned.fastq.gz",
        ref=config["RNA"]["ITS"]["Fungi_UNITE_db"]
    output:
        its_r1=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/ITS_candidates/{{sample}}_ITS_candidates_1.fastq.gz",
        its_r2=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/ITS_candidates/{{sample}}_ITS_candidates_2.fastq.gz",
    log:
        stdout=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/logs/bbduk_its.log"
    benchmark:
        f"{RNA_INTERMEDIATE_DIR}/{{sample}}/benchmarks/bbduk_its.txt"
    params:
        options=config["RNA"]["ITS"]["bbduk"]["options"]
    threads:
        config["RNA"]["ITS"]["bbduk"]["threads"]
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.its_r1})

        bbduk.sh \
            in={input.cleaned_r1} \
            in2={input.cleaned_r2} \
            outm={output.its_r1} \
            outm2={output.its_r2} \
            ref={input.ref} \
            threads={threads} \
            {params.options} \
            > {log.stdout} 2>&1
        """