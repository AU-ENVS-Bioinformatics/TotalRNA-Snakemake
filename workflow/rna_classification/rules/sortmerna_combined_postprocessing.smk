from pathlib import Path


RNA_ENVS_DIR = str(Path(workflow.basedir) / "RNA" / "envs")


if (
    config["RNA"]["method"] == "sortmerna"
    and config["RNA"]["refinement"] == "combined"
):

    # =====================================================
    # Extract aligned, SSU and LSU read IDs from SAM
    # =====================================================

    rule extract_sortmerna_combined_ids:
        conda:
            f"{RNA_ENVS_DIR}/seqkit.yaml"
        message:
            "[SortMeRNA] Extracting rRNA, SSU and LSU read IDs for {wildcards.sample}"
        input:
            sam=rules.sortmerna_combined.output.sam
        output:
            alignment_info=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/sortmerna/combined/{{sample}}_alignment_info.tsv",
            rrna_ids=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/sortmerna/combined/{{sample}}_rRNA_ids.txt",
            ssu_ids=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/sortmerna/combined/{{sample}}_SSU_ids.txt"
        log:
            stdout=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/logs/sortmerna_combined_ids.log"
        benchmark:
            f"{RNA_INTERMEDIATE_DIR}/{{sample}}/benchmarks/sortmerna_combined_ids.txt"
        shell:
            r"""
            set -euo pipefail

            mkdir -p "$(dirname "{output.alignment_info}")"
            mkdir -p "$(dirname "{log.stdout}")"

            samtools view -F 4 "{input.sam}" \
                | awk -v OFS="\t" \
                    '{{print $1, $3, $6, $12, $13, length($10)}}' \
                > "{output.alignment_info}" \
                2> "{log.stdout}"

            cut -f1 "{output.alignment_info}" \
                | sort -u \
                > "{output.rrna_ids}"

            awk '$2 ~ /^SSU_/' "{output.alignment_info}" \
                | cut -f1 \
                | sort -u \
                > "{output.ssu_ids}"
            """

    # =====================================================
    # Extract non-rRNA reads
    # =====================================================

    rule extract_non_rRNA_sortmerna_combined:
        conda:
            f"{RNA_ENVS_DIR}/seqkit.yaml"
        message:
            "[SortMeRNA] Extracting non-rRNA reads "
            "for {wildcards.sample}"
        input:
            cleaned_r1=f"{QC_DIR}/{{sample}}/decontamination/{{sample}}_R1.cleaned.fastq.gz",
            cleaned_r2=f"{QC_DIR}/{{sample}}/decontamination/{{sample}}_R2.cleaned.fastq.gz",
            read_ids=rules.extract_sortmerna_combined_ids.output.rrna_ids
        output:
            r1=f"{RNA_CLASSIFIED_DIR}/{{sample}}/non_rRNA/{{sample}}_non_rRNA_1.fastq.gz",
            r2=f"{RNA_CLASSIFIED_DIR}/{{sample}}/non_rRNA/{{sample}}_non_rRNA_2.fastq.gz"
        log:
            stdout=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/logs/extract_sortmerna_combined_non_rRNA.log"
        benchmark:
            f"{RNA_INTERMEDIATE_DIR}/{{sample}}/benchmarks/extract_sortmerna_combined_non_rRNA.txt"
        threads:
            config["RNA"]["threads"]
        shell:
            r"""
            set -euo pipefail

            mkdir -p "$(dirname "{output.r1}")"
            mkdir -p "$(dirname "{log.stdout}")"

            seqkit grep \
                -v \
                -f "{input.read_ids}" \
                --threads {threads} \
                "{input.cleaned_r1}" \
                -o "{output.r1}" \
                2> "{log.stdout}"

            seqkit grep \
                -v \
                -f "{input.read_ids}" \
                --threads {threads} \
                "{input.cleaned_r2}" \
                -o "{output.r2}" \
                2>> "{log.stdout}"
            """
    
    # =====================================================
    # Extract SSU reads
    # =====================================================

    rule extract_SSU_sortmerna_combined:
        conda:
            f"{RNA_ENVS_DIR}/seqkit.yaml"
        message:
            "[SortMeRNA] Extracting competitive SSU reads for {wildcards.sample}"
        input:
            cleaned_r1=f"{QC_DIR}/{{sample}}/decontamination/{{sample}}_R1.cleaned.fastq.gz",
            cleaned_r2=f"{QC_DIR}/{{sample}}/decontamination/{{sample}}_R2.cleaned.fastq.gz",
            read_ids=rules.extract_sortmerna_combined_ids.output.ssu_ids
        output:
            r1=f"{RNA_CLASSIFIED_DIR}/{{sample}}/SSU/{{sample}}_SSU_1.fastq.gz",
            r2=f"{RNA_CLASSIFIED_DIR}/{{sample}}/SSU/{{sample}}_SSU_2.fastq.gz"
        log:
            stdout=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/logs/extract_sortmerna_combined_SSU.log"
        benchmark:
            f"{RNA_INTERMEDIATE_DIR}/{{sample}}/benchmarks/extract_sortmerna_combined_SSU.txt"
        threads:
            config["RNA"]["threads"]
        shell:
            r"""
            set -euo pipefail

            mkdir -p "$(dirname "{output.r1}")"
            mkdir -p "$(dirname "{log.stdout}")"

            seqkit grep \
                -f "{input.read_ids}" \
                --threads {threads} \
                "{input.cleaned_r1}" \
                -o "{output.r1}" \
                2> "{log.stdout}"

            seqkit grep \
                -f "{input.read_ids}" \
                --threads {threads} \
                "{input.cleaned_r2}" \
                -o "{output.r2}" \
                2>> "{log.stdout}"
            """