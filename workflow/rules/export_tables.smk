rule export_tsv:
    input:
        db="results/{folder}/database.db",
        reads_done="results/{folder}/read_length.done",
    output:
        "results/{folder}/tsv/mapped_reads_{model}_{threshold}.tsv",
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/{folder}/export_tsv_{model}_{threshold}.log",
    benchmark:
        "results/benchmarks/export_tsv_{folder}_{model}_{threshold}.txt",
    params:
        model=lambda wc: wc.model,
        threshold=lambda wc: wc.threshold,
    script:
        "../scripts/export_tsv.py"