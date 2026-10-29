rule bracken:
    conda:
        "../envs/bracken.yaml"
    wildcard_constraints:
        marker="SSU|ITS"
    message:
        "[Bracken] Re-estimating {wildcards.marker} abundance for {wildcards.sample} against {wildcards.database}"
    input:
        kraken_report=f"{TAXONOMY_DIR}/{{sample}}/{{marker}}/{{database}}/kraken/{{sample}}.{{database}}.k2report"
    output:
        bracken_output=f"{TAXONOMY_DIR}/{{sample}}/{{marker}}/{{database}}/bracken/{{sample}}.{{database}}.bracken.tsv",
        bracken_kreport_output=f"{TAXONOMY_DIR}/{{sample}}/{{marker}}/{{database}}/bracken/{{sample}}.{{database}}.bracken.k2report",
    log:
        stdout=f"{TAXONOMY_DIR}/{{sample}}/{{marker}}/{{database}}/logs/{{sample}}.{{database}}.bracken.log",
    benchmark:
        f"{TAXONOMY_DIR}/{{sample}}/{{marker}}/{{database}}/benchmarks/{{sample}}.{{database}}.bracken.txt",
    params:
        database=lambda wc: kraken_db(wc, config=config),
        level=lambda wc: bracken_level(wc, config=config),
        options=lambda wc: bracken_options(wc, config=config),
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bracken_output}) $(dirname {log.stdout})
        skip=0
        if [ ! -s {input.kraken_report} ]; then
            echo "[Bracken] Empty Kraken report; no abundance estimates." > {log.stdout}
            skip=1
        else
            status=0
            bracken \
                -d {params.database} \
                -i {input.kraken_report} \
                -o {output.bracken_output} \
                -w {output.bracken_kreport_output} \
                {params.level} \
                {params.options} \
                > {log.stdout} 2>&1 || status=$?
            if [ "$status" -ne 0 ]; then
                if grep -Fxq "Error: no reads found. Please check your Kraken report" {log.stdout}; then
                    echo "[WARNING] Bracken produced no estimates at the requested rank. Original Kraken classifications are retained separately." >> {log.stdout}
                    skip=1
                else
                    exit "$status"
                fi
            fi
        fi
        if [ "$skip" -eq 1 ]; then
            printf 'name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads\n' > {output.bracken_output}
            : > {output.bracken_kreport_output}
        elif [ ! -s {output.bracken_output} ] || [ ! -s {output.bracken_kreport_output} ]; then
            echo "[ERROR] Bracken exited successfully but expected outputs are missing or empty." >> {log.stdout}
            exit 1
        fi
        """

