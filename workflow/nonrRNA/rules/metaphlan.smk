rule metaphlan:
    conda:
        "/data_2/rasmus/conda-envs/humann_env"
    message:
        "[metaphlan] species-level microbial profiling for {wildcards.sample}"
    input:
        nonrRNA_concatenate=rules.concatenate.output.nonrRNA_concatenate
    output:
        metaphlan_profile=f"{RESULTS_DIR}/nonrRNA/{{sample}}/metaphlan/{{sample}}_metaphlan_profile.tsv",
    log:
        stdout=f"{RESULTS_DIR}/nonrRNA/{{sample}}/logs/metaphlan.log"
    benchmark:
        f"{RESULTS_DIR}/nonrRNA/{{sample}}/benchmarks/metaphlan.txt"
    params:
        db=config["databases"]["humann_bowtie2db"],
        options=config["nonrRNA"]["metaphlan"]["options"]
    threads:
        config["nonrRNA"]["metaphlan"].get("threads", 8)
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.metaphlan_profile})
    
        metaphlan {input.nonrRNA_concatenate} \
            {params.options} \
            --nproc {threads} \
            --bowtie2db {params.db} \
            --output_file {output.metaphlan_profile} \
            > {log.stdout} 2>&1
        """
#metaphlan ANN11_concantenated.fastq.gz --input_type fastq --bowtie2db /data_2/Databases/humann/metaphlan/ --nproc 16 --index mpa_vJun23_CHOCOPhlAnSGB_202307 --output_file ../humann/Metaphlan_Jun23_profile.txt