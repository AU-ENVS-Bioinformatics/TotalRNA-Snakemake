
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


    rule link_LSU_sortmerna:
        message:
            "[SortMeRNA] Publishing LSU reads "
            "for {wildcards.sample}"

        input:
            r1=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/sortmerna/{FINAL_STAGE_NAME}/{{sample}}_{FINAL_STAGE_NAME}.aligned_fwd.fq.gz",
            r2=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/sortmerna/{FINAL_STAGE_NAME}/{{sample}}_{FINAL_STAGE_NAME}.aligned_rev.fq.gz",

        output:
            r1=f"{RNA_CLASSIFIED_DIR}/{{sample}}/LSU/{{sample}}_LSU_1.fastq.gz",
            r2=f"{RNA_CLASSIFIED_DIR}/{{sample}}/LSU/{{sample}}_LSU_2.fastq.gz",

        shell:
            r"""
            set -euo pipefail

            mkdir -p "$(dirname "{output.r1}")"

            ln -sfn "{input.r1}" "{output.r1}"
            ln -sfn "{input.r2}" "{output.r2}"
            """


    rule link_non_SSU_sortmerna:
        message:
            "[SortMeRNA] Publishing non-SSU rRNA reads "
            "for {wildcards.sample}"

        input:
            r1=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/sortmerna/{FINAL_STAGE_NAME}/{{sample}}_{FINAL_STAGE_NAME}.aligned_fwd.fq.gz",
            r2=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/sortmerna/{FINAL_STAGE_NAME}/{{sample}}_{FINAL_STAGE_NAME}.aligned_rev.fq.gz",

        output:
            r1=f"{RNA_CLASSIFIED_DIR}/{{sample}}/non_SSU/{{sample}}_non_SSU_1.fastq.gz",
            r2=f"{RNA_CLASSIFIED_DIR}/{{sample}}/non_SSU/{{sample}}_non_SSU_2.fastq.gz",

        shell:
            r"""
            set -euo pipefail

            mkdir -p "$(dirname "{output.r1}")"

            ln -sfn "{input.r1}" "{output.r1}"
            ln -sfn "{input.r2}" "{output.r2}"
            """


    rule combine_rRNA_sortmerna:
        message:
            "[SortMeRNA] Combining SSU and LSU reads "
            "for {wildcards.sample}"

        input:
            ssu_r1=rules.link_SSU_sortmerna.output.r1,
            ssu_r2=rules.link_SSU_sortmerna.output.r2,
            lsu_r1=rules.link_LSU_sortmerna.output.r1,
            lsu_r2=rules.link_LSU_sortmerna.output.r2,

        output:
            r1=f"{RNA_CLASSIFIED_DIR}/{{sample}}/rRNA/{{sample}}_rRNA_1.fastq.gz",
            r2=f"{RNA_CLASSIFIED_DIR}/{{sample}}/rRNA/{{sample}}_rRNA_2.fastq.gz",

        shell:
            r"""
            set -euo pipefail

            mkdir -p "$(dirname "{output.r1}")"

            cat "{input.ssu_r1}" "{input.lsu_r1}" > "{output.r1}"
            cat "{input.ssu_r2}" "{input.lsu_r2}" > "{output.r2}"
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