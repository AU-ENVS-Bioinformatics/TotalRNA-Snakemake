from pathlib import Path

ENVS_DIR = Path(workflow.basedir) / "taxonomy/envs"

rule soil_kraken_biom:
    conda:
        f"{ENVS_DIR}/kraken_biom.yaml"
    message:
        "[Kraken-Biom] Combining Soil Bracken reports against {wildcards.database}"
    input:
        reports=habitat_taxonomy_reports_for_database
    output:
        biom=f"{HABITAT_ANALYSIS_DIR}/summary/Soil/{{database}}/taxonomy/taxonomy.biom"
    log:
        stdout=f"{HABITAT_ANALYSIS_DIR}/summary/Soil/{{database}}/taxonomy/logs/kraken_biom.log"
    benchmark:
        f"{HABITAT_ANALYSIS_DIR}/summary/Soil/{{database}}/taxonomy/benchmarks/kraken_biom.txt"
    params:
        options=config["Habitat_refinement"]["Soil"]["kraken_biom"]["options"],
        level=lambda wc: config["Habitat_refinement"]["Soil"]["kraken_biom"][wc.database]["level"],
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.biom})

        kraken-biom {input.reports} \
            -o {output.biom} \
            {params.level} \
            {params.options}
        """
    
rule soil_kraken_biom_convert:
    conda:
        f"{ENVS_DIR}/kraken_biom.yaml"
    message:
        "[biom convert] converting {wildcards.database} BIOM to TSV"
    input:
        biom=f"{HABITAT_ANALYSIS_DIR}/summary/Soil/{{database}}/taxonomy/taxonomy.biom"
    output:
        biom_tsv=f"{HABITAT_ANALYSIS_DIR}/summary/Soil/{{database}}/taxonomy/taxonomy.tsv"
    shell:
        r"""
        mkdir -p $(dirname {output.biom_tsv})

        biom convert -i {input.biom} -o {output.biom_tsv} --to-tsv --header-key taxonomy
        """

rule soil_biom_to_phyloseq_raw:
    conda:
        f"{ENVS_DIR}/r_phyloseq.yaml"
    message:
        "[phyloseq] Creating raw phyloseq object from {wildcards.database}"
    input:
        biom=f"{HABITAT_ANALYSIS_DIR}/summary/Soil/{{database}}/taxonomy/taxonomy.biom"
    output:
        rds=f"{HABITAT_ANALYSIS_DIR}/summary/Soil/{{database}}/phyloseq/phyloseq.rds"
    shell:
        r"""
        mkdir -p $(dirname {output.rds})

        Rscript -e '
        suppressPackageStartupMessages({{
            library(biomformat)
            library(phyloseq)
        }})

        physeq <- import_biom("{input.biom}")

        saveRDS(
            physeq,
            file="{output.rds}",
            compress="xz"
        )
        '
        """