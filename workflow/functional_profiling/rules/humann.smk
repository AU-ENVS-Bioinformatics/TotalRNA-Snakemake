rule concat_sample_pe_reads:
    conda:
        "../envs/seqtk.yaml"
    message:
        "[Concatenate] Concatenate paired-end reads for {wildcards.sample}"
    input:
        nonrna_r1=f"{RNA_CLASSIFIED_DIR}/{{sample}}/non_rRNA/{{sample}}_non_rRNA_1.fastq.gz",
        nonrna_r2=f"{RNA_CLASSIFIED_DIR}/{{sample}}/non_rRNA/{{sample}}_non_rRNA_2.fastq.gz",
    output:
        pe_merged_reads=f"{RNA_CLASSIFIED_DIR}/{{sample}}/non_rRNA/{{sample}}_non_rRNA_pe_merged.fastq.gz",
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.pe_merged_reads})

        seqtk mergepe {input.nonrna_r1} {input.nonrna_r2} | gzip > {output.pe_merged_reads}
        """

rule humann:
    conda:
        "../envs/humann.yaml"
    message:
        "[humann] profiling the presence/absence and abundance of microbial pathways in a community for {wildcards.sample}"
    input:
        merged_pe_reads=rules.concat_sample_pe_reads.output.pe_merged_reads,
        metaphlan_profile=rules.metaphlan.output.metaphlan_profile
    output:
        genefamilies=f"{FUNCTION_DIR_READS}/{{sample}}/humann/{{sample}}_genefamilies.tsv",
        abundance=f"{FUNCTION_DIR_READS}/{{sample}}/humann/{{sample}}_pathabundance.tsv",
        coverage=f"{FUNCTION_DIR_READS}/{{sample}}/humann/{{sample}}_pathcoverage.tsv"
    log:
        stdout=f"{FUNCTION_DIR_READS}/{{sample}}/logs/humann.log",
        time=f"{FUNCTION_DIR_READS}/{{sample}}/benchmarks/humann_time.txt"
    benchmark:
        f"{FUNCTION_DIR_READS}/{{sample}}/benchmarks/humann.txt"
    params:
        result_dir=f"{FUNCTION_DIR_READS}/{{sample}}/humann/",
        nucleotide_db=config["functional_profiling"]["humann"]["nucleotide_database"],
        protein_db=config["functional_profiling"]["humann"]["protein_database"],
        basename=f"{{sample}}",
        options=config["functional_profiling"]["humann"]["options"]
    threads:
        config["functional_profiling"]["humann"].get("threads", 8)
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.genefamilies})

        /usr/bin/time -v -o {log.time} \
            humann --input {input.merged_pe_reads} \
                --output {params.result_dir} \
                --threads {threads} \
                --output-basename {params.basename} \
                --taxonomic-profile {input.metaphlan_profile} \
                --nucleotide-database {params.nucleotide_db} \
                --protein-database {params.protein_db} \
                {params.options} \
                > {log.stdout} 2>&1
        """
