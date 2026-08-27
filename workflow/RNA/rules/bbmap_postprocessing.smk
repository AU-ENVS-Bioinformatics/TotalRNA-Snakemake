
rule SSU_rRNA_ids:
    conda:
        "../envs/bbmap.yaml"
    message:
        "[BBMap] identify aligned rRNA reads for {wildcards.sample}"
    input:
        all_bam_reads=rules.bbmap.output.all_bam_reads,
    output:
        rna_seqid=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/bbmap/{{sample}}_rRNA_readid.txt",
    log:
        stdout=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/logs/bbmap_seqread_ids.log"
    benchmark:
        f"{RNA_INTERMEDIATE_DIR}/{{sample}}/benchmarks/bbmap_seqread_ids.txt"
    params:
        filters="-F 4 -q 20", # consider adding additional filters to ensure we get high confidence rRNA reads, such as higher mapping quality, or specific flags to ensure we only get primary alignments, and not secondary or supplementary alignments.
        prefix="SSU" # filter for specific rRNA types based on reference database, e.g. SSU, LSU, 5S, etc. This can be adjusted based on the specific reference database used and the types of rRNA you want to focus on.
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.rna_seqid})

        samtools view {params.filters} {input.all_bam_reads} \
            | grep '{params.prefix}' \
            | cut -f1 | cut -d' ' -f1 | sort | uniq \
            > {output.rna_seqid} 2> {log.stdout}
        """


rule non_rRNA_ids:
    conda:
        "../envs/bbmap.yaml"
    message:
        "[BBMap] identify aligned non-rRNA reads for {wildcards.sample}"
    input:
        all_bam_reads=rules.bbmap.output.all_bam_reads,
    output:
        non_rna_seqid=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/bbmap/{{sample}}_non_rRNA_readid.txt",
    log:
        stdout=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/logs/bbmap_seqread_ids.log"

    benchmark:
        f"{RNA_INTERMEDIATE_DIR}/{{sample}}/benchmarks/bbmap_seqread_ids.txt"
    output:
        non_rna_seqid=f"{RESULTS_DIR}/RNA/{{sample}}/bbmap/{{sample}}_non_rRNA_readid.txt",
    log:
        stdout=f"{RESULTS_DIR}/RNA/{{sample}}/logs/bbmap_seqread_ids.log"
    benchmark:
        f"{RESULTS_DIR}/RNA/{{sample}}/benchmarks/bbmap_seqread_ids.log"
    params:
        filters="-F 4"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.non_rna_seqid})

        samtools view {params.filters} {input.all_bam_reads} \
            | cut -f1 | cut -d' ' -f1 | sort | uniq > {output.non_rna_seqid} 2> {log.stdout}
        """

rule link_rRNA_bbmap:
    conda:
        "../envs/bbmap.yaml"
    message:
        "[BBMap] filter aligned rRNA reads for {wildcards.sample}"
    input:
        cleaned_r1=f"{QC_DIR}/{{sample}}/decontamination/{{sample}}_R1.cleaned.fastq.gz",
        cleaned_r2=f"{QC_DIR}/{{sample}}/decontamination/{{sample}}_R2.cleaned.fastq.gz",
        read_ids=rules.SSU_rRNA_ids.output.rna_seqid
    output:
        filtered_r1=f"{RNA_CLASSIFIED_DIR}/{{sample}}/rRNA/{{sample}}_rRNA_1.fastq.gz",
        filtered_r2=f"{RNA_CLASSIFIED_DIR}/{{sample}}/rRNA/{{sample}}_rRNA_2.fastq.gz",
    log:
        stdout=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/logs/bbmap_rRNA_extraction.log"
    benchmark:
        f"{RNA_INTERMEDIATE_DIR}/{{sample}}/benchmarks/bbmap_rRNA_extraction.txt"
    threads:
        config["RNA"]["threads"]
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.filtered_r1})

        seqkit grep \
            -f {input.read_ids} \
            --threads {threads} \
            {input.cleaned_r1} \
            -o {output.filtered_r1} \
            2> {log.stdout}

        seqkit grep \
            -f {input.read_ids} \
            --threads {threads} \
            {input.cleaned_r2} \
            -o {output.filtered_r2} \
            2>> {log.stdout}
        """

rule link_non_rRNA_bbmap:
    conda:
        "../envs/bbmap.yaml"
    message:
        "[BBMap] filter aligned non-rRNA reads for {wildcards.sample}"
    input:
        cleaned_r1=f"{QC_DIR}/{{sample}}/decontamination/{{sample}}_R1.cleaned.fastq.gz",
        cleaned_r2=f"{QC_DIR}/{{sample}}/decontamination/{{sample}}_R2.cleaned.fastq.gz",
        read_ids=rules.non_rRNA_ids.output.non_rna_seqid
    output:
        filtered_r1=f"{RNA_CLASSIFIED_DIR}/{{sample}}/non_rRNA/{{sample}}_non_rRNA_1.fastq.gz",
        filtered_r2=f"{RNA_CLASSIFIED_DIR}/{{sample}}/non_rRNA/{{sample}}_non_rRNA_2.fastq.gz",
    log:
        stdout=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/logs/bbmap_non_rRNA_extraction.log"
    benchmark:
        f"{RNA_INTERMEDIATE_DIR}/{{sample}}/benchmarks/bbmap_non_rRNA_extraction.txt"
    threads:
        config["RNA"]["threads"]
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.filtered_r1})

        seqkit grep \
            -f {input.read_ids} \
            --threads {threads} \
            {input.cleaned_r1} \
            -o {output.filtered_r1} \
            2> {log.stdout}

        seqkit grep \
            -f {input.read_ids} \
            --threads {threads} \
            {input.cleaned_r2} \
            -o {output.filtered_r2} \
            2>> {log.stdout}
        """