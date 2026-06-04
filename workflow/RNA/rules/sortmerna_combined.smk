
#consider changing the mkdir -p and workdir to have a basename = dirname {output.sam} and then reuse it, such that we dont have to have same "dir" in both params and outpiút
rule sortmerna_combined:
    conda:
        "../envs/sortmerna.yaml"
    message:
        "[SortMeRNA] seperate reads  for {wildcards.sample} according to the numerous databases at once"
    input:
        cleaned_r1=f"{RESULTS_DIR}/qc/{{sample}}/decontamination/{{sample}}_R1.cleaned.fastq.gz",
        cleaned_r2=f"{RESULTS_DIR}/qc/{{sample}}/decontamination/{{sample}}_R2.cleaned.fastq.gz",
        db=config["databases"][f"sortmeRNA_ssu_lsu"],
        db_idx=config["databases"][f"sortmeRNA_ssu_lsu_idx"]
    output:
        sam=f"{RESULTS_DIR}/RNA/{{sample}}/sortmerna/aligned/{{sample}}_ssu_lsu.aligned.sam",
    log:
        stdout = f"{RESULTS_DIR}/RNA/{{sample}}/logs/sortmerna_aligned.log"
    benchmark:
        f"{RESULTS_DIR}/RNA/{{sample}}/benchmarks/sortmerna_aligned.txt"
    params:
        workdir=f"{RESULTS_DIR}/RNA/{{sample}}/sortmerna/aligned",
        options=config["RNA"]["sortmerna_combined"]["options"],
        aligned_prefix=f"{RESULTS_DIR}/RNA/{{sample}}/sortmerna/aligned/{{sample}}_ssu_lsu.aligned",
        nonaligned_prefix=f"{RESULTS_DIR}/RNA/{{sample}}/sortmerna/aligned/{{sample}}_ssu_lsu.nonaligned"
    threads:
        config["RNA"]["sortmerna_combined"]["threads"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.workdir}

        sortmerna \
            --ref {input.db} \
            --idx-dir {input.db_idx} \
            --workdir {params.workdir} \
            --threads {threads} \
            --reads {input.cleaned_r1} \
            --reads {input.cleaned_r2} \
            {params.options} \
            --aligned {params.aligned_prefix} \
            > {log.stdout} 2>&1
        
        rm -rf {params.workdir}/kvdb || true
        rm -rf {params.workdir}/readb || true
        """

#sortmerna --ref /data_2/Databases/SILVA_138/SILVA_138.1_LSU_SSU_Ref_NR99_tax_silva_trunc.fasta --idx-dir /data_2/Databases/sortmerna_idx/SSU_LSU/idx/ 
#--workdir sortmerna/ --reads decontamination/ANN_11_R1.cleaned.fastq.gz --reads decontamination/ANN_11_R2.cleaned.fastq.gz 
#--sam --SQ --print_all_reads -zip-out 0 --threads 12
