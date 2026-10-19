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

        mkdir -p $(dirname {output.bracken_output})

        bracken \
          -d {params.database} \
          -i {input.kraken_report} \
          -o {output.bracken_output} \
          -w {output.bracken_kreport_output} \
          {params.level} \
          {params.options} \
          > {log.stdout} 2>&1
        """