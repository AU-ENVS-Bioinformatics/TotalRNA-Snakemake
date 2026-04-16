rule fastqc:
    conda:
        "../envs/fastqc.yaml"
    message:
        "[FastQC] {wildcards.sample} {wildcards.readtag}"
    input:
        fq=lambda wc: fastqc_input_path(wc.sample, wc.readtag)
    output:
        qcdir=directory(
            "{output_dir}/{sample}/QC/fastqc/{readtag}"
        )
    log:
        stdout = "{output_dir}/{sample}/logs/fastqc_{readtag}.log"
    benchmark:
        "{output_dir}/{sample}/benchmarks/fastqc_{readtag}.txt"
    threads:
        config["qc"]["fastqc"].get("threads", 2)
    wildcard_constraints:
        readtag="R1|R2|SE"
    shell:
        r"""
        mkdir -p {output.qcdir} \
                 "$(dirname {log.stdout})"

        fastqc \
            --threads {threads} \
            --outdir {output.qcdir} \
            "{input.fq}" \
            > {log.stdout} 2>&1
        """
