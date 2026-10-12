rule metaphlan:
    conda:
        "../envs/humann.yaml"
    message:
        "[metaphlan] species-level microbial profiling for {wildcards.sample}"
    input:
        nonrna_r1=f"{RNA_CLASSIFIED_DIR}/{{sample}}/non_rRNA/{{sample}}_non_rRNA_1.fastq.gz",
        nonrna_r2=f"{RNA_CLASSIFIED_DIR}/{{sample}}/non_rRNA/{{sample}}_non_rRNA_2.fastq.gz",
    output:
        metaphlan_profile=f"{FUNCTION_DIR_READS}/{{sample}}/metaphlan/{{sample}}_metaphlan_profile.tsv",
    log:
        stdout=f"{FUNCTION_DIR_READS}/{{sample}}/logs/metaphlan.log",
        time=f"{FUNCTION_DIR_READS}/{{sample}}/benchmarks/metaphlan_time.txt"
    benchmark:
        f"{FUNCTION_DIR_READS}/{{sample}}/benchmarks/metaphlan.txt"
    params:
        db=config["databases"]["humann_bowtie2db"],
        options=config["functional_profiling"]["metaphlan"]["options"],
        metaphlan_bowtie2=f"{FUNCTION_DIR_READS}/{{sample}}/metaphlan/{{sample}}_metaphlan_profile.tsv.bowtie2.bz2",
    threads:
        config["functional_profiling"]["metaphlan"].get("threads", 8)
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.metaphlan_profile})

        rm -f {params.metaphlan_bowtie2}

        /usr/bin/time -v -o {log.time} \
            metaphlan "{input.nonrna_r1},{input.nonrna_r2}" \
                --nproc {threads} \
                --bowtie2db {params.db} \
                --bowtie2out {params.metaphlan_bowtie2} \
                --output_file {output.metaphlan_profile} \
                {params.options} \
                > {log.stdout} 2>&1
        """
