
rule SSU_rRNA_ids:
    conda:
        "../envs/bbmap.yaml"
    message:
        "[BBMap] identify aligned ssu reads for {wildcards.sample}"
    input:
        all_RNA_reads_bam=rules.bbmap.output.all_RNA_reads_bam,
    output:
        ssu_seqid=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/bbmap/SSU/{{sample}}_ssu_readid.txt",
    log:
        stdout=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/logs/bbmap_ssu_read_ids.log"
    benchmark:
        f"{RNA_INTERMEDIATE_DIR}/{{sample}}/benchmarks/bbmap_ssu_read_ids.txt"
    params:
        prefix="SSU", # filter for specific rRNA types based on reference database, e.g. SSU, LSU, 5S, etc. This can be adjusted based on the specific reference database used and the types of rRNA you want to focus on.
        filters="-F 4" # -f 3 1 read paired (0x1) + 2  read mapped in proper pair (0x2) but this exclude potential SSU reads, where the mates align to different contigs, which due to conserved regions might be a plausible scenario. Therefore, we will use -F 4 to exclude unmapped reads and match on prefix.
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.ssu_seqid})
        mkdir -p $(dirname {log.stdout})

        samtools view {params.filters} {input.all_RNA_reads_bam} \
            | grep '{params.prefix}' \
            | cut -f1 | cut -d' ' -f1 | sort | uniq \
            > {output.ssu_seqid} 2> {log.stdout}
        """


rule non_rRNA_ids:
    conda:
        "../envs/bbmap.yaml"
    message:
        "[BBMap] identifying nonaligned SSU and LSU reads resulting in non-rRNA reads for {wildcards.sample}"
    input:
        all_RNA_reads_bam=rules.bbmap.output.all_RNA_reads_bam,
    output:
        non_rna_seqid=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/bbmap/non_rRNA/{{sample}}_non_rRNA_readid.txt",
    log:
        stdout=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/logs/bbmap_nonrRNA_ids.log"
    benchmark:
        f"{RNA_INTERMEDIATE_DIR}/{{sample}}/benchmarks/bbmap_nonrRNA_ids.txt"
    output:
        non_rna_seqid=f"{RESULTS_DIR}/RNA/{{sample}}/bbmap/non_rRNA/{{sample}}_non_rRNA_readid.txt",
    log:
        stdout=f"{RESULTS_DIR}/RNA/{{sample}}/logs/bbmap_nonrRNA_ids.log"
    benchmark:
        f"{RESULTS_DIR}/RNA/{{sample}}/benchmarks/bbmap_nonrRNA_ids.txt"
    params:
        filters="-f 13" #1 read paired (0x1) + 4 read unmapped (0x4) + 8 mate unmapped (0x8)
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.non_rna_seqid})

        samtools view {params.filters} {input.all_RNA_reads_bam} \
            | cut -f1 | cut -d' ' -f1 | sort | uniq > {output.non_rna_seqid} 2> {log.stdout}
        """

rule bbmap_ssu:
    conda:
        "../envs/bbmap.yaml"
    message:
        "[BBmap] Extracting SSU reads for {wildcards.sample}"
    input:
        cleaned_r1=f"{QC_DIR}/{{sample}}/decontamination/{{sample}}_R1.cleaned.fastq.gz",
        cleaned_r2=f"{QC_DIR}/{{sample}}/decontamination/{{sample}}_R2.cleaned.fastq.gz",
        read_ids=rules.SSU_rRNA_ids.output.ssu_seqid
    output:
        filtered_r1=f"{RNA_CLASSIFIED_DIR}/{{sample}}/SSU/{{sample}}_SSU_1.fastq.gz",
        filtered_r2=f"{RNA_CLASSIFIED_DIR}/{{sample}}/SSU/{{sample}}_SSU_2.fastq.gz",
    log:
        stdout=f"{RNA_CLASSIFIED_DIR}/{{sample}}/logs/bbmap_ssu_extraction.log"
    benchmark:
        f"{RNA_CLASSIFIED_DIR}/{{sample}}/benchmarks/bbmap_ssu_extraction.txt"
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

rule bbmap_non_rrna:
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
        stdout=f"{RNA_CLASSIFIED_DIR}/{{sample}}/logs/bbmap_non_rRNA_extraction.log"
    benchmark:
        f"{RNA_CLASSIFIED_DIR}/{{sample}}/benchmarks/bbmap_non_rRNA_extraction.txt"
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