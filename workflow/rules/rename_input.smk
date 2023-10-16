rule rename:
    output:
        touch("results/rename.done"),
    log:
        "logs/renaming_files.log",
    conda:
        "../envs/base_python.yaml"
    script:
        "../scripts/rename_files.py"


rule skip_rename:
    output:
        touch("results/rename.done"),
    log:
        "logs/renaming_files.log",
    conda:
        "../envs/base_python.yaml"
    script:
        "../scripts/skip_rename_files.py"
