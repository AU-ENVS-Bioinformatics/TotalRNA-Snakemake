rule fastqc:
    conda:
        "../envs/fastqc.yaml"
    message:
        "[FastQC] quality check for {wildcards.sample} {wildcards.read}"
    input:
        fq=lambda wc: raw_r1(wc) if wc.read == "R1" else raw_r2(wc),
    output:
        directory(
            f"{QC_DIR}/{{sample}}/fastqc/{{read}}"
        )
    log:
        stdout=f"{QC_DIR}/{{sample}}/logs/fastqc_{{read}}.log"
    benchmark:
        f"{QC_DIR}/{{sample}}/benchmarks/fastqc_{{read}}.txt"
    threads:
        config["qc"]["fastqc"].get("threads", 2)
    wildcard_constraints:
        read="R1|R2"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {output}

        fastqc \
          -t {threads} \
          -o {output} \
          {input.fq} \
          > {log.stdout} 2>&1
        """