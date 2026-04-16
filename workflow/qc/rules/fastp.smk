rule fastp:
    conda:
        "../envs/fastp.yaml"
    log:
        stdout=lambda wc: f"{logdir(wc)}/fastp.log"
    benchmark:
        lambda wc: f"{benchmarkdir(wc)}/fastp.txt"
    message:
        "[Fastp]: Running fastp on {wildcards.sample}"
    input:
        r1=lambda wc: fastqs(wc)[0],
        r2=lambda wc: fastqs(wc)[1]
    output:
        trim_r1="results/preprocessing/{sample}/fastp/{sample}_R1.trim.fastq.gz",
        trim_r2="results/preprocessing/{sample}/fastp/{sample}_R2.trim.fastq.gz",
        trim_json="results/preprocessing/{sample}/fastp/{sample}_fastp.json",
        trim_html="results/preprocessing/{sample}/fastp/{sample}_fastp.html",
    params:
        options = lambda wc: config["fastp"].get("options", ""),
    threads:
        config["fastp"].get("threads", 8)
    shell:
        """
        outdir=$(dirname {output.trim_r1})
        mkdir -p "$outdir"

        fastp \
          --in1 {input.r1} --in2 {input.r2} \
          --out1 {output.trim_r1} --out2 {output.trim_r2} \
            --html {output.trim_html} --json {output.trim_json} \
                --thread {threads} {params.options} \
                     > {log.stdout} 2>&1

        """