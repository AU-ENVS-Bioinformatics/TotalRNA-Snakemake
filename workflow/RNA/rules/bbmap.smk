rule bbmap:
    conda:
        "../envs/bbmap.yaml"
    message:
        "[BBMap] Running BBMap for {wildcards.sample}"
    input:
        cleaned_r1=f"{RESULTS_DIR}/qc/{{sample}}/decontamination/{{sample}}_R1.cleaned.fastq.gz",
        cleaned_r2=f"{RESULTS_DIR}/qc/{{sample}}/decontamination/{{sample}}_R2.cleaned.fastq.gz",
        ref=config["databases"]["sortmeRNA_ssu_lsu"]
    output:
        all_bam_reads=f"{RESULTS_DIR}/RNA/{{sample}}/bbmap/{{sample}}_all_reads.bam",
    log:
        stdout = f"{RESULTS_DIR}/RNA/{{sample}}/logs/bbmap.log"
    benchmark:
        f"{RESULTS_DIR}/RNA/{{sample}}/benchmarks/bbmap.txt"
    params:
        options=config["RNA"]["bbmap"]["options"]
    threads:
        config["RNA"]["bbmap"]["threads"]
    shell:
        r"""
        set -euo pipefail
        
        bbmap.sh \
            in={input.cleaned_r1} \
            in2={input.cleaned_r2} \
            ref={input.ref} \
            out={output.all_bam_reads} \
            threads={threads} \
            {params.options} \
            > {log.stdout} 2>&1
        """
#ref=/data_2/Databases/SILVA_138/SILVA_138.1_LSU_SSU_Ref_NR99_tax_silva_trunc.fasta in=ANN_11/decontamination/ANN_11_R1.cleaned.fastq.gz in2=ANN_11/decontamination/ANN_11_R2.cleaned.fastq.gz minid=0.90 ambiguous=best unpigz=t maxindel=100 minhits=2 threads=16 out=test2.bam
