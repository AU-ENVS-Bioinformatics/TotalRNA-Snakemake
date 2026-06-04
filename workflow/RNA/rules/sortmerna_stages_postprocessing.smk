
if config["RNA"]["method"] == "sortmerna" and config["RNA"]["refinement"] == "staged":
    rule link_rRNA_sortmerna:
        message:
            "[SortMeRNA] linking rRNA from first stage for {wildcards.sample}"
        input:
            r1=lambda wc: f"{RESULTS_DIR}/RNA/{wc.sample}/sortmerna/{FIRST_STAGE_NAME}/{wc.sample}_{FIRST_STAGE_NAME}.aligned_fwd.fq.gz",
            r2=lambda wc: f"{RESULTS_DIR}/RNA/{wc.sample}/sortmerna/{FIRST_STAGE_NAME}/{wc.sample}_{FIRST_STAGE_NAME}.aligned_rev.fq.gz",
        output:
            r1=f"{RESULTS_DIR}/rRNA/{{sample}}/filtered/{{sample}}_SSU.rRNA.R1.fastq.gz",
            r2=f"{RESULTS_DIR}/rRNA/{{sample}}/filtered/{{sample}}_SSU.rRNA.R2.fastq.gz",
        shell:
            """
            mkdir -p $(dirname {output.r1})
            ln -sf {input.r1} {output.r1}
            ln -sf {input.r2} {output.r2}
            """

    rule link_non_rRNA_sortmerna:
        message:
            "[SortMeRNA] linking final non-rRNA for {wildcards.sample}"
        input:
            r1=lambda wc: f"{RESULTS_DIR}/RNA/{wc.sample}/sortmerna/{FINAL_STAGE_NAME}/{wc.sample}_{FINAL_STAGE_NAME}.nonaligned_fwd.fq.gz",
            r2=lambda wc: f"{RESULTS_DIR}/RNA/{wc.sample}/sortmerna/{FINAL_STAGE_NAME}/{wc.sample}_{FINAL_STAGE_NAME}.nonaligned_rev.fq.gz",
        output:
            r1=f"{RESULTS_DIR}/nonrRNA/{{sample}}/filtered/{{sample}}_R1.nonrRNA.fastq.gz",
            r2=f"{RESULTS_DIR}/nonrRNA/{{sample}}/filtered/{{sample}}_R2.nonrRNA.fastq.gz",
        shell:
            """
            mkdir -p $(dirname {output.r1})
            ln -sf {input.r1} {output.r1}
            ln -sf {input.r2} {output.r2}
            """

#sortmerna --ref /data_2/Databases/SILVA_138/SILVA_138.1_LSU_SSU_Ref_NR99_tax_silva_trunc.fasta --idx-dir /data_2/Databases/sortmerna_idx/SSU_LSU --workdir sortmerna/ --reads decontamination/ANN_11_R1.cleaned.fastq.gz --reads decontamination/ANN_11_R2.cleaned.fastq.gz --sam --SQ --print_all_reads -zip-out 0 --threads 12