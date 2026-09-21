from pathlib import Path

HABITAT_ENVS_DIR = Path(workflow.basedir) / "taxonomy/envs"

rule habitat_kraken_biom:
    conda:
        f"{HABITAT_ENVS_DIR}/kraken_biom.yaml"
    message:
        "[Kraken-Biom] Combining Bracken reports for {wildcards.habitat}/{wildcards.database}"
    input:
        reports=habitat_reports_for_database
    output:
        biom=f"{HABITAT_ANALYSIS_DIR}/summary/{{habitat}}/{{database}}/taxonomy/taxonomy.biom"
    log:
        stdout=f"{HABITAT_ANALYSIS_DIR}/summary/{{habitat}}/{{database}}/taxonomy/logs/kraken_biom.log"
    benchmark:
        f"{HABITAT_ANALYSIS_DIR}/summary/{{habitat}}/{{database}}/taxonomy/benchmarks/kraken_biom.txt"
    params:
        options=habitat_biom_options,
        level=habitat_biom_level,
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.biom})

        kraken-biom {input.reports} \
          -o {output.biom} \
          {params.level} \
          {params.options} \
          > {log.stdout} 2>&1
        """


rule habitat_kraken_biom_convert:
    conda:
        f"{HABITAT_ENVS_DIR}/kraken_biom.yaml"
    message:
        "[biom convert] Converting {wildcards.habitat}/{wildcards.database} BIOM to TSV"
    input:
        biom=f"{HABITAT_ANALYSIS_DIR}/summary/{{habitat}}/{{database}}/taxonomy/taxonomy.biom"
    output:
        biom_tsv=f"{HABITAT_ANALYSIS_DIR}/summary/{{habitat}}/{{database}}/taxonomy/taxonomy.tsv"
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.biom_tsv})

        biom convert \
          -i {input.biom} \
          -o {output.biom_tsv} \
          --to-tsv \
          --header-key taxonomy
        """


rule habitat_biom_to_phyloseq_raw:
    conda:
        f"{HABITAT_ENVS_DIR}/r_phyloseq.yaml"
    message:
        "[phyloseq] Creating raw phyloseq object for {wildcards.habitat}/{wildcards.database}"
    input:
        biom=f"{HABITAT_ANALYSIS_DIR}/summary/{{habitat}}/{{database}}/taxonomy/taxonomy.biom"
    output:
        rds=f"{HABITAT_ANALYSIS_DIR}/summary/{{habitat}}/{{database}}/phyloseq/phyloseq.rds"
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.rds})

        Rscript -e '
        suppressPackageStartupMessages({{
            library(biomformat)
            library(phyloseq)
        }})

        physeq <- import_biom("{input.biom}")

        saveRDS(
            physeq,
            file = "{output.rds}",
            compress = "xz"
        )
        '
        """