rule fastp:
    input:
        R1=lambda wc: sample_to_illumina[wc.sample][0],
        R2=lambda wc: sample_to_illumina[wc.sample][1]
    output:
        trim_R1=f"{illumina_trim_read_path}/{{sample}}_1.trim.fastq.gz",
        trim_R2=f"{illumina_trim_read_path}/{{sample}}_2.trim.fastq.gz",
        trim_html=f"{illumina_trim_read_path}/{{sample}}_trim.html",
        trim_json=f"{illumina_trim_read_path}/{{sample}}_trim.json"
    log:
        "../logs/{sample}/fastp.log"
    benchmark:
        "../benchmarks/{sample}/fastp.tsv"
    conda:
        "../envs/fastp.yaml"
    message:
        "[Fastp] Trimming reads for sample {wildcards.sample}"
    threads:
        config["fastp"].get("threads", 8)
    resources:
        mem_mb=config["fastp"].get("mem_mb", 4000),
    params:
        options=config["fastp"].get("options", "")
    shell:
        """
        outdir=$(dirname {output.trim_R1})
        mkdir -p "$outdir"
        
        fastp \
            --in1 {input.R1} \
            --in2 {input.R2} \
            --out1 {output.trim_R1} \
            --out2 {output.trim_R2} \
            --html {output.trim_html} \
            --json {output.trim_json} \
            --thread {threads} \
            {params.options} \
            &> {log}
        """