rule habitat_kraken2:
    conda:
        "../envs/kraken2.yaml"
    message:
        "[Kraken2] Classifying {wildcards.habitat} reads from {wildcards.sample} against {wildcards.database}"
    input:
        read_1=habitat_read_1,
        read_2=habitat_read_2,
    output:
        report=f"{HABITAT_ANALYSIS_DIR}/{{sample}}/{{habitat}}/{{database}}/kraken/{{sample}}.{{database}}.k2report",
        kraken=f"{HABITAT_ANALYSIS_DIR}/{{sample}}/{{habitat}}/{{database}}/kraken/{{sample}}.{{database}}.kraken",
    log:
        stdout=f"{HABITAT_ANALYSIS_DIR}/{{sample}}/{{habitat}}/{{database}}/logs/{{sample}}.{{database}}.kraken2.log",
    benchmark:
        f"{HABITAT_ANALYSIS_DIR}/{{sample}}/{{habitat}}/{{database}}/benchmarks/{{sample}}.{{database}}.kraken2.txt"
    params:
        db=habitat_kraken_db,
        options=habitat_kraken_options,
    threads:
        habitat_kraken_threads
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.report})

        kraken2 \
          --db {params.db} \
          --threads {threads} \
          --report {output.report} \
          --output {output.kraken} \
          --paired \
          {input.read_1} \
          {input.read_2} \
          {params.options} \
          > "{log.stdout}" 2>&1
        """