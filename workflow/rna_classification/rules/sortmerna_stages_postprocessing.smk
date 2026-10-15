
if config["RNA"]["method"] == "sortmerna" and config["RNA"]["refinement"] == "staged":
    rule link_SSU_sortmerna:
        message:
            "[SortMeRNA] Publishing SSU reads "
            "for {wildcards.sample}"
        input:
            r1=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/sortmerna/{FIRST_STAGE_NAME}/{{sample}}_{FIRST_STAGE_NAME}.aligned_fwd.fq.gz",
            r2=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/sortmerna/{FIRST_STAGE_NAME}/{{sample}}_{FIRST_STAGE_NAME}.aligned_rev.fq.gz",
        output:
            r1=f"{RNA_CLASSIFIED_DIR}/{{sample}}/SSU/{{sample}}_SSU_1.fastq.gz",
            r2=f"{RNA_CLASSIFIED_DIR}/{{sample}}/SSU/{{sample}}_SSU_2.fastq.gz",
        shell:
            r"""
            set -euo pipefail

            mkdir -p "$(dirname "{output.r1}")"

            ln -sfn "{input.r1}" "{output.r1}"
            ln -sfn "{input.r2}" "{output.r2}"
            """

    rule link_non_rRNA_sortmerna:
        message:
            "[SortMeRNA] Publishing final non-rRNA reads "
            "for {wildcards.sample}"
        input:
            r1=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/sortmerna/{FINAL_STAGE_NAME}/{{sample}}_{FINAL_STAGE_NAME}.nonaligned_fwd.fq.gz",
            r2=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/sortmerna/{FINAL_STAGE_NAME}/{{sample}}_{FINAL_STAGE_NAME}.nonaligned_rev.fq.gz",
        output:
            r1=f"{RNA_CLASSIFIED_DIR}/{{sample}}/non_rRNA/{{sample}}_non_rRNA_1.fastq.gz",
            r2=f"{RNA_CLASSIFIED_DIR}/{{sample}}/non_rRNA/{{sample}}_non_rRNA_2.fastq.gz",
        shell:
            r"""
            set -euo pipefail

            mkdir -p "$(dirname "{output.r1}")"

            ln -sfn "{input.r1}" "{output.r1}"
            ln -sfn "{input.r2}" "{output.r2}"
            """


# we dont really need to link LSU reads, since we are not using them for anything, but we can do it if we want to have them for downstream analysis. For now, we will just link the SSU and non-rRNA reads, since those are the ones we are interested in.

#    rule link_LSU_sortmerna:
#        message:
#            "[SortMeRNA] Publishing LSU reads for {wildcards.sample}"
#        input:
#            r1=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/sortmerna/{FINAL_STAGE_NAME}/{{sample}}_{FINAL_STAGE_NAME}.aligned_fwd.fq.gz",
#            r2=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/sortmerna/{FINAL_STAGE_NAME}/{{sample}}_{FINAL_STAGE_NAME}.aligned_rev.fq.gz",
#        output:
#            r1=f"{RNA_CLASSIFIED_DIR}/{{sample}}/LSU/{{sample}}_LSU_1.fastq.gz",
#            r2=f"{RNA_CLASSIFIED_DIR}/{{sample}}/LSU/{{sample}}_LSU_2.fastq.gz",
#        shell:
#            r"""
#            set -euo pipefail
#
#            mkdir -p "$(dirname "{output.r1}")"
#
#            ln -sfn "{input.r1}" "{output.r1}"
#            ln -sfn "{input.r2}" "{output.r2}"
#            """
