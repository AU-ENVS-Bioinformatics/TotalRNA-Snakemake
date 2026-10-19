rule habitat_kraken2:
    conda:
        "../envs/kraken2.yaml"
    message:
        "[Kraken2] Classifying {wildcards.habitat} reads from {wildcards.sample} against {wildcards.database}"
    input:
        read_1=lambda wc: habitat_read_1(
            wc,
            results_dir=RESULTS_DIR,
        ),
        read_2=lambda wc: habitat_read_2(
            wc,
            results_dir=RESULTS_DIR,
        ),
    output:
        report=f"{HABITAT_ANALYSIS_DIR}/{{sample}}/{{habitat}}/{{database}}/kraken/{{sample}}.{{database}}.k2report",
        kraken=f"{HABITAT_ANALYSIS_DIR}/{{sample}}/{{habitat}}/{{database}}/kraken/{{sample}}.{{database}}.kraken",
    log:
        stdout=f"{HABITAT_ANALYSIS_DIR}/{{sample}}/{{habitat}}/{{database}}/logs/{{sample}}.{{database}}.kraken2.log",
    benchmark:
        f"{HABITAT_ANALYSIS_DIR}/{{sample}}/{{habitat}}/{{database}}/benchmarks/{{sample}}.{{database}}.kraken2.txt"
    params:
        db=lambda wc: habitat_kraken_db(
            wc,
            config=config,
        ),
        options=lambda wc: habitat_kraken_options(
            wc,
            config=config,
        ),
    threads:
        lambda wc: habitat_kraken_threads(
            wc,
            config=config,
        )
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