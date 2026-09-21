rule habitat_bracken:
    conda:
        "../envs/bracken.yaml"
    message:
        "[Bracken] Re-estimating abundance for {wildcards.sample} using {wildcards.habitat}/{wildcards.database}"
    input:
        kraken_report=f"{HABITAT_ANALYSIS_DIR}/{{sample}}/{{habitat}}/{{database}}/kraken/{{sample}}.{{database}}.k2report"
    output:
        bracken_output=f"{HABITAT_ANALYSIS_DIR}/{{sample}}/{{habitat}}/{{database}}/bracken/{{sample}}.{{database}}.bracken.tsv",
        bracken_kreport_output=f"{HABITAT_ANALYSIS_DIR}/{{sample}}/{{habitat}}/{{database}}/bracken/{{sample}}.{{database}}.bracken.k2report",
    log:
        stdout=f"{HABITAT_ANALYSIS_DIR}/{{sample}}/{{habitat}}/{{database}}/logs/{{sample}}.{{database}}.bracken.log",
    benchmark:
            f"{HABITAT_ANALYSIS_DIR}/{{sample}}/{{habitat}}/{{database}}/benchmarks/{{sample}}.{{database}}.bracken.txt"
    params:
        db=habitat_kraken_db,
        level=habitat_bracken_level,
        options=habitat_bracken_options,
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.bracken_output})

        bracken \
          -d {params.db} \
          -i {input.kraken_report} \
          -o {output.bracken_output} \
          -w {output.bracken_kreport_output} \
          {params.level} \
          {params.options} \
          > {log.stdout} 2>&1
        """