trim_galore_params = config.get("trim_galore", "")
trim_galore_threads = config["threads"]["trim_galore"]
trim_galore_params.append(f"--cores {trim_galore_threads}")


rule trim_galore_pe:
    input:
        r1=lambda wc: samples_dict[wc.sample].r1,
        r2=lambda wc: samples_dict[wc.sample].r2,
    output:
        fasta_fwd="results/trim_galore/{sample}_fwd.fq.gz",
        fasta_rev="results/trim_galore/{sample}_rev.fq.gz",
    threads: trim_galore_threads
    params:
        extra=" ".join(trim_galore_params),
        outdir="results/trim_galore",
        report_dir="results/trim_galore/reports"
    log:
        "logs/trim_galore/{sample}.log"
    benchmark:
        "results/benchmarks/trim_galore_{sample}.txt"
    conda:
        "../envs/trim_galore.yaml"
    shell:
        r"""
        mkdir -p {params.outdir} {params.report_dir}

        trim_galore \
            --paired \
            --basename {wildcards.sample} \
            {params.extra} \
            --output_dir {params.outdir} \
            {input.r1} {input.r2} \
            > {log} 2>&1

        mv {params.outdir}/{wildcards.sample}*_val_1.fq.gz {output.fasta_fwd} >> {log} 2>&1
        mv {params.outdir}/{wildcards.sample}*_val_2.fq.gz {output.fasta_rev} >> {log} 2>&1
        mv {params.outdir}/{wildcards.sample}*_trimming_report.txt {params.report_dir} >> {log} 2>&1
        """