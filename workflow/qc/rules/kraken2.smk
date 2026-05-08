rule kraken2:
    conda:
        "../envs/kraken2.yaml"
    message:
        "[Kraken2] estimate species composition for {wildcards.sample}"
    input:
        trim_r1="{outdir}/{sample}/QC/trimmed/{sample}_R1.fastq.gz",
        trim_r2="{outdir}/{sample}/QC/trimmed/{sample}_R2.fastq.gz",
    output:
        report="{outdir}/{sample}/QC/classification/{sample}.report",
        kraken="{outdir}/{sample}/QC/classification/{sample}.kraken"
    log:
        stdout = "{outdir}/{sample}/logs/kraken2.log"
    benchmark:
        "{outdir}/{sample}/benchmarks/kraken2.txt"
    params:
        db=config["databases"]["kraken_db"],
        options=config["qc"]["kraken2"]["options"],
    threads:
        config["qc"]["kraken2"].get("threads", 2)
    wildcard_constraints:
        outdir=".+"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.report})

        kraken2 \
          --db {params.db} \
          --threads {threads} \
          --report {output.report} \
          --output {output.kraken} \
          {params.options} \
          {input.trim_r1} {input.trim_r2} \
          > {log.stdout} 2>&1
        """

rule kraken2_parse:
    conda:
        "../../envs/python.yaml"
    message:
        "[Parse Kraken Report] Contamination check for {wildcards.sample}"
    input:
        report="{outdir}/{sample}/QC/classification/{sample}.report"
    output:
        flag="{outdir}/{sample}/QC/classification/{sample}.human.flag"
    log:
        stdout = "{outdir}/{sample}/logs/kraken2_parse.log"
    benchmark:
        "{outdir}/{sample}/benchmarks/kraken2_parse.txt"
    params:
        threshold=config["qc"]["kraken2_parse"]["human_threshold"],
        taxid=config["qc"]["kraken2_parse"]["taxid"],
    wildcard_constraints:
        outdir=".+"
    run:
        from qc.scripts.qc_utils import (
            parse_kraken_human_percentage,
            decide_flag,
        )

        pct = parse_kraken_human_percentage(input.report, params.taxid)
        decision = decide_flag(pct, params.threshold)

        with open(output.flag, "w") as f:
            f.write(decision + "\n")

        with open(log.stdout, "w") as l:
            l.write(f"Human percentage: {pct}\nDecision: {decision}\n")

#https://github.com/jenniferlu717/KrakenTools/blob/master/extract_kraken_reads.py
rule remove_human_reads:
    conda:
        "../envs/kraken2.yaml"
    message:
        "[Exclude kraken2 reads] extract and exclude contaimation reads based on taxid for {wildcards.sample}"
    input:
        kraken="{outdir}/{sample}/QC/classification/{sample}.kraken",
        report="{outdir}/{sample}/QC/classification/{sample}.report",
        flag="{outdir}/{sample}/QC/classification/{sample}.human.flag",
        trim_r1="{outdir}/{sample}/QC/trimmed/{sample}_R1.fastq.gz",
        trim_r2="{outdir}/{sample}/QC/trimmed/{sample}_R2.fastq.gz",
    output:
        cleaned_r1="{outdir}/{sample}/QC/decontamination/{sample}_R1.clean.fastq.gz",
        cleaned_r2="{outdir}/{sample}/QC/decontamination/{sample}_R2.clean.fastq.gz",
    log:
        stdout = "{outdir}/{sample}/logs/kraken2_exclude.log"
    benchmark:
        "{outdir}/{sample}/benchmarks/kraken2_exclude.txt"
    params:
        threshold=config["qc"]["kraken2_parse"]["human_threshold"],
        taxid=config["qc"]["kraken2_parse"]["taxid"],
    wildcard_constraints:
        outdir=".+"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.cleaned_r1})

        decision=$(cat {input.flag})

        echo "Decision: $decision"

        if [ "$decision" = "REMOVE" ]; then
            echo "Removing reads"
            extract_kraken_reads.py \
                -k {input.kraken} \
                --report {input.report} \
                -s1 {input.trim_r1} \
                -s2 {input.trim_r2} \
                -t {params.taxid} \
                --include-children \
                --exclude \
                --fastq-output \
                -o {output.cleaned_r1} \
                -o2 {output.cleaned_r2} 
        else
            echo "Keeping original reads"
            cp {input.trim_r1} {output.cleaned_r1}
            cp {input.trim_r2} {output.cleaned_r2}
        fi
        """