rule export_tsv:
    input:
        db="results/{folder}/database.db",
        reads_done="results/{folder}/read_length.done",
    output:
        "results/{folder}/tsv/mapped_reads_exported.tsv",
        "results/{folder}/tsv/read_length_exported.tsv",
        "results/{folder}/tsv/contig_length_exported.tsv"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/{folder}/export_tsv.log",
    benchmark:
        "results/benchmarks/export_tsv_{folder}.txt",
    script:
        "../scripts/export_tsv.py"


