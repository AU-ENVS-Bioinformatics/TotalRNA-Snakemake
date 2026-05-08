rule bwa_decontamination:
    conda:
        "../../envs/bwa.yaml"
    message:
        "[BWA] decontaminate {wildcards.sample}"
    log:
        stdout = "{out}/{sample}/logs/bwa_decontamination.log"
    benchmark:
        "{out}/{sample}/benchmarks/bwa_decontamination.txt"
    input:
        r1=trimmed_r1,
        r2=trimmed_r2,
        flag="{out}/{sample}/QC/classification/{sample}.human.flag"
    output:
        r1=clean_r1,
        r2=clean_r2
    params:
        genome=config["databases"]["human_genome"],
        out=lambda wc: f"{outdir(wc.sample)}/{wc.sample}/QC/decontaminated"
    threads: config["qc"]["bwa_decontamination"]["threads"]
    run:
        from workflow.utils.qc_logic import load_decontam_flag

        remove = load_decontam_flag(input.flag)

        shell(f"mkdir -p {params.out} $(dirname {log})")

        if remove:
            shell(f"""
                bwa mem -t {threads} {params.genome} \
                {input.r1} {input.r2} \
                | samtools view -b -f 4 \
                | samtools fastq \
                -1 {output.r1} \
                -2 {output.r2} \
                > {log} 2>&1
            """)
        else:
            shell(f"""
                ln -sf {input.r1} {output.r1}
                ln -sf {input.r2} {output.r2}
                echo "No contamination detected" > {log}
            """)