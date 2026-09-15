rule kraken2_soil:
    conda:
        "../envs/kraken2.yaml"
    message:
        "[Kraken2:] Classifying paired non-rRNA reads from {wildcards.sample} against the habitat-specific database {wildcards.database}"
    input:
        non_rna_r1=f"{RNA_CLASSIFIED_DIR}/{{sample}}/non_rRNA/{{sample}}_non_rRNA_2.fastq.gz",
        non_rna_r2=f"{RNA_CLASSIFIED_DIR}/{{sample}}/non_rRNA/{{sample}}_non_rRNA_2.fastq.gz",
    output:
        report=f"{HABITAT_ANALYSIS_DIR}/{{sample}}/Soil/{{database}}/kraken/{{sample}}.{{database}}.k2report",
        kraken=f"{HABITAT_ANALYSIS_DIR}/{{sample}}/Soil/{{database}}/kraken/{{sample}}.{{database}}.kraken"
    log:
        stdout=f"{HABITAT_ANALYSIS_DIR}/{{sample}}/Soil/{{database}}/logs/{{sample}}.{{database}}.kraken2.log"
    benchmark:
        f"{HABITAT_ANALYSIS_DIR}/{{sample}}/Soil/{{database}}/benchmarks/{{sample}}.{{database}}.kraken2.txt"
    params:
        db=lambda wildcards: config["Habitat_refinement"]["Soil"]["kraken2"]["databases"][wildcards.database]["path"],
        options=config["Habitat_refinement"]["Soil"]["kraken2"].get("options","")
    threads:
        config["Habitat_refinement"]["Soil"]["kraken2"].get("threads",2)
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.report})

        kraken2 \
          --db {params.db} \
          --threads {threads} \
          --report {output.report} \
          --output {output.kraken} \
          --paired {input.non_rna_r1} {input.non_rna_r2} \
            {params.options} \
          > {log.stdout} 2>&1
        """
