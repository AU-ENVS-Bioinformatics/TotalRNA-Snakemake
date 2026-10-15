
#consider changing the mkdir -p and workdir to have a basename = dirname {output.sam} and then reuse it, such that we dont have to have same "dir" in both params and outpiút
rule sortmerna_combined:
    conda:
        "../envs/sortmerna.yaml"
    message:
        "[SortMeRNA] seperate reads  for {wildcards.sample} according to the numerous databases at once"
    input:
        cleaned_r1=f"{QC_DIR}/{{sample}}/decontamination/{{sample}}_R1.cleaned.fastq.gz",
        cleaned_r2=f"{QC_DIR}/{{sample}}/decontamination/{{sample}}_R2.cleaned.fastq.gz",
        db=config["databases"][f"sortmeRNA_ssu_lsu"],
        db_idx=config["databases"][f"sortmeRNA_ssu_lsu_idx"]
    output:
        sam=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/sortmerna/combined/{{sample}}_ssu_lsu.aligned.sam"
    log:
        stdout=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/logs/sortmerna_combined.log"
    benchmark:
        f"{RNA_INTERMEDIATE_DIR}/{{sample}}/benchmarks/sortmerna_combined.txt"
    params:
        workdir=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/sortmerna/combined",
        aligned_prefix=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/sortmerna/combined/{{sample}}_ssu_lsu.aligned",
        options=config["RNA"]["sortmerna_combined"]["options"]
    threads:
        config["RNA"]["sortmerna_combined"]["threads"]
    shell:
        r"""
        set -euo pipefail

        mkdir -p "{params.workdir}"
        mkdir -p "$(dirname "{log.stdout}")"

        sortmerna \
            --ref "{input.db}" \
            --idx-dir "{input.db_idx}" \
            --workdir "{params.workdir}" \
            --threads {threads} \
            --reads "{input.cleaned_r1}" \
            --reads "{input.cleaned_r2}" \
            --aligned "{params.aligned_prefix}" \
            {params.options} \
            > "{log.stdout}" 2>&1

        rm -rf {params.workdir}/kvdb || true
        rm -rf {params.workdir}/readb || true
        """

#sortmerna --ref /data_2/Databases/SILVA_138/SILVA_138.1_LSU_SSU_Ref_NR99_tax_silva_trunc.fasta --idx-dir /data_2/Databases/sortmerna_idx/SSU_LSU/idx/ 
#--workdir sortmerna/ --reads decontamination/ANN_11_R1.cleaned.fastq.gz --reads decontamination/ANN_11_R2.cleaned.fastq.gz 
#--sam --SQ --print_all_reads -zip-out 0 --threads 12
