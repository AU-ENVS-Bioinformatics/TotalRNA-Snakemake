from pathlib import Path

ENVS_DIR = Path(workflow.basedir) / "taxonomy/envs"
SCRIPTS_DIR = Path(workflow.basedir) / "taxonomy/scripts"

rule kraken_biom:
    conda:
        f"{ENVS_DIR}/kraken_biom.yaml"
    message:
        "[Kraken-Biom] Combining SSU reports across samples against {wildcards.database}"
    input:
        reports=taxonomy_reports_for_database
    output:
        biom=f"{TAXONOMY_DIR}/summary/SSU/{{database}}/taxonomy/taxonomy.biom"
    log:
        stdout=f"{TAXONOMY_DIR}/summary/SSU/{{database}}/logs/kraken_biom.log"
    benchmark:
        f"{TAXONOMY_DIR}/summary/SSU/{{database}}/benchmarks/kraken_biom.txt"
    params:
        options=config["taxonomy"]["SSU"]["kraken_biom"]["options"],
        level=biom_level,
    shell:
        r"""
        mkdir -p $(dirname {output.biom})

        kraken-biom {input.reports} \
            -o {output.biom} \
            {params.level} \
            {params.options} \
            > {log.stdout} 2>&1
        """

#kraken-biom *.k2report -o bracken_genus.biom --min S --max D
rule biom_table_stats_qualitative:
    conda:
        f"{ENVS_DIR}/kraken_biom.yaml"
    message:
        "[biom summarize-table] Creating qualitative SSU summary for against {wildcards.database}"
    input:
        biom=f"{TAXONOMY_DIR}/summary/SSU/{{database}}/taxonomy/taxonomy.biom"
    output:
        biom_qual=f"{TAXONOMY_DIR}/summary/SSU/{{database}}/statistics/taxonomy_qualitative.txt"
    log:
        stdout=f"{TAXONOMY_DIR}/summary/SSU/{{database}}/logs/biom_summary_qualitative.log"
    benchmark:
        f"{TAXONOMY_DIR}/summary/SSU/{{database}}/benchmarks/biom_summary_qualitative.txt"
    shell:
        r"""
        mkdir -p $(dirname {output.biom_qual})

        biom summarize-table -i {input.biom} \
            -o {output.biom_qual} \
            --qualitative \
            > {log.stdout} 2>&1
        """

rule biom_table_stats_observations:
    conda:
        f"{ENVS_DIR}/kraken_biom.yaml"
    message:
        "[biom summarize-table] Creating SSU observation summary for {wildcards.database}"
    input:
        biom=f"{TAXONOMY_DIR}/summary/SSU/{{database}}/taxonomy/taxonomy.biom"
    output:
        biom_obs=f"{TAXONOMY_DIR}/summary/SSU/{{database}}/statistics/taxonomy_observations.txt"
    log:
        stdout=f"{TAXONOMY_DIR}/summary/SSU/{{database}}/logs/biom_summary_observations.log"
    benchmark:
        f"{TAXONOMY_DIR}/summary/SSU/{{database}}/benchmarks/biom_summary_observations.txt"
    shell:
        r"""
        mkdir -p $(dirname {output.biom_obs})

        biom summarize-table -i {input.biom} \
            -o {output.biom_obs} \
            --observations \
            > {log.stdout} 2>&1
        """

rule kraken_biom_convert:
    conda:
        f"{ENVS_DIR}/kraken_biom.yaml"
    message:
        "[biom convert] converting {wildcards.database} SSU BIOM to TSV"
    input:
        biom=f"{TAXONOMY_DIR}/summary/SSU/{{database}}/taxonomy/taxonomy.biom"
    output:
        biom_tsv=f"{TAXONOMY_DIR}/summary/SSU/{{database}}/taxonomy/taxonomy.tsv"
    log:
        f"{TAXONOMY_DIR}/summary/SSU/{{database}}/logs/biom_convert.log"
    benchmark:
        f"{TAXONOMY_DIR}/summary/SSU/{{database}}/benchmarks/biom_convert.txt"
    shell:
        r"""
        mkdir -p $(dirname {output.biom_tsv})

        biom convert -i {input.biom} -o {output.biom_tsv} --to-tsv --header-key taxonomy > {log} 2>&1
        """

rule biom_to_phyloseq_raw:
    conda:
        f"{ENVS_DIR}/r_phyloseq.yaml"
    message:
        "[phyloseq] Creating raw SSU phyloseq object from {wildcards.database}"
    input:
        biom=f"{TAXONOMY_DIR}/summary/SSU/{{database}}/taxonomy/taxonomy.biom"
    output:
        rds=f"{TAXONOMY_DIR}/summary/SSU/{{database}}/phyloseq/phyloseq_raw.rds"
    log:
        f"{TAXONOMY_DIR}/summary/SSU/{{database}}/logs/biom_to_phyloseq_raw.log"
    benchmark:
        f"{TAXONOMY_DIR}/summary/SSU/{{database}}/benchmarks/biom_to_phyloseq_raw.txt"
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
        ' > {log} 2>&1
        """

rule biom_to_phyloseq_filtered:
    conda:
        f"{ENVS_DIR}/r_phyloseq.yaml"
    message:
        "[phyloseq] Creating filtered SSU phyloseq object from {wildcards.database}, by removing multi-cellular organism "
    input:
        biom=f"{TAXONOMY_DIR}/summary/SSU/{{database}}/taxonomy/taxonomy.biom"
    output:
        rds=f"{TAXONOMY_DIR}/summary/SSU/{{database}}/phyloseq/phyloseq_filtered.rds"
    log:
        f"{TAXONOMY_DIR}/summary/SSU/{{database}}/logs/biom_to_phyloseq.log"
    benchmark:
        f"{TAXONOMY_DIR}/summary/SSU/{{database}}/benchmarks/biom_to_phyloseq.txt"
    params:
        prefix=f"{TAXONOMY_DIR}/summary/SSU/{{database}}/"
    shell:
        r"""
        Rscript {SCRIPTS_DIR}/biom_to_phyloseq.R {input.biom} {output.rds} {params.prefix} > {log} 2>&1
        """

#mamba create -n kraken_biom_test -c conda-forge -c bioconda python=3.10 kraken-biom biom-format