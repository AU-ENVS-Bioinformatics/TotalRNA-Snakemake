rule rename:
    output:
        touch("results/rename.done"),
    log:
        "logs/renaming_files.log",
    conda:
        "../envs/base_python.yaml"
    benchmark:
        "benchmarks/renaming_files.log"
    script:
        "../scripts/rename_files.py"
