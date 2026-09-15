rule kraken2:
    conda:
        "../envs/kraken2.yaml"
    wildcard_constraints:
        marker="SSU|ITS"
    message:
        "[Kraken2] Classifying {wildcards.marker} reads from {wildcards.sample} against {wildcards.database}"
    input:
        read_1=taxonomy_read_1,
        read_2=taxonomy_read_2,
    output:
        report=f"{TAXONOMY_DIR}/{{sample}}/{{marker}}/{{database}}/kraken/{{sample}}.{{database}}.k2report",
        kraken=f"{TAXONOMY_DIR}/{{sample}}/{{marker}}/{{database}}/kraken/{{sample}}.{{database}}.kraken",
    log:
        stdout=f"{TAXONOMY_DIR}/{{sample}}/{{marker}}/{{database}}/logs/{{sample}}.{{database}}.kraken2.log",
    benchmark:
        f"{TAXONOMY_DIR}/{{sample}}/{{marker}}/{{database}}/benchmarks/{{sample}}.{{database}}.kraken2.txt",
    params:
        db=kraken_db,
        options=kraken_options,
    threads:
        kraken_threads
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.kraken})

        kraken2 \
          --db {params.db} \
          --threads {threads} \
          --report {output.report} \
          --output {output.kraken} \
          --paired \
          {input.read_1} \
          {input.read_2} \
          {params.options} \
          > {log.stdout} 2>&1
        """