rule bracken_soil:
    conda:
        "../envs/bracken.yaml"
    message:
        "[Bracken] re-estimate abundance for {wildcards.sample} for {wildcards.database}"
    input:
        kraken_report=f"{HABITAT_ANALYSIS_DIR}/{{sample}}/Soil/{{database}}/kraken/{{sample}}.{{database}}.k2report",
    output:
        bracken_output=f"{HABITAT_ANALYSIS_DIR}/{{sample}}/Soil/{{database}}/bracken/{{sample}}.{{database}}.bracken.tsv",
        bracken_kreport_output=f"{HABITAT_ANALYSIS_DIR}/{{sample}}/Soil/{{database}}/bracken/{{sample}}.{{database}}.bracken.k2report"
    log:
        stdout=f"{HABITAT_ANALYSIS_DIR}/{{sample}}/Soil/{{database}}/logs/{{sample}}.{{database}}.bracken.log"
    benchmark:
        f"{HABITAT_ANALYSIS_DIR}/{{sample}}/Soil/{{database}}/benchmarks/{{sample}}.{{database}}.bracken.txt"
    params:
        db=lambda wildcards: config["Habitat_refinement"]["Soil"]["kraken2"]["databases"][wildcards.database]["path"],
        level=lambda wildcards: config["Habitat_refinement"]["Soil"]["bracken"][wildcards.database]["level"],
        options=config["Habitat_refinement"]["Soil"]["bracken"]["options"]
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