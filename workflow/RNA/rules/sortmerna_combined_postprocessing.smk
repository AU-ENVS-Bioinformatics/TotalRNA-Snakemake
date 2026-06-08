from pathlib import Path

ENVS_DIR = Path(workflow.basedir)/"RNA/envs"
SCRIPTS_DIR = Path(workflow.basedir)/"RNA/scripts"

if config["RNA"]["method"] == "sortmerna" and config["RNA"]["refinement"] == "combined":
    
    rule extract_all_aligned_ids:
        input:
            sam=rules.sortmerna_combined.output.sam
        output:
            readid_info_tsv=f"{RESULTS_DIR}/RNA/{{sample}}/sortmerna/aligned/{{sample}}.all_aligned_ids.tsv",
            readid_txt=f"{RESULTS_DIR}/RNA/{{sample}}/sortmerna/aligned/{{sample}}.all_aligned_ids.txt",
        shell:  
            r"""
            set -euo pipefail

            # full TSV for downstream filtering
            samtools view {input.sam} | awk -v OFS="\t" '{{print $1,$3,$6,$12,$13,length($10)}}' | tee {output.readid_info_tsv} | cut -f1 | sort -u > {output.readid_txt}
            """
    
    rule filter_and_assign_taxonomy:
        conda:
            f"{ENVS_DIR}/rna_python_tools.yaml"
        message:
                "[Sortmerna_postprocessing_filter] for {wildcards.sample} filtering reads of percentage of the read length and at least certain number of nucleotide matches to the reference, and extract the ENA accession number for each read id and assign taxonomy to each read based on the ncbi taxonomy table"
        input:
            id_info=rules.extract_all_aligned_ids.output.readid_info_tsv,
            taxonomy=config["databases"]["taxonomy_table"]
        output:
            tsv=f"{RESULTS_DIR}/RNA/{{sample}}/sortmerna/aligned/{{sample}}_SSU_filtered_taxonomy.tsv",
            read_ids=f"{RESULTS_DIR}/RNA/{{sample}}/sortmerna/aligned/{{sample}}_SSU_filtered_taxonomy_ids.txt",
            failed_ids=f"{RESULTS_DIR}/RNA/{{sample}}/sortmerna/aligned/{{sample}}_failed_alignment_ids.txt",
            discarded_ids=f"{RESULTS_DIR}/RNA/{{sample}}/sortmerna/aligned/{{sample}}_discarded_alignment_ids.txt",
            taxonomy_discarded_id=f"{RESULTS_DIR}/RNA/{{sample}}/sortmerna/aligned/{{sample}}_discarded_taxonomy_ids.txt",
            statistics=f"{RESULTS_DIR}/RNA/{{sample}}/sortmerna/aligned/{{sample}}_statistics.tsv"
        log:
            stdout = f"{RESULTS_DIR}/RNA/{{sample}}/logs/sortmerna_postprocessing_filter_taxonomy.log"
        benchmark:
            f"{RESULTS_DIR}/RNA/{{sample}}/benchmarks/sortmerna_postprocessing_filter_taxonomy.txt"
        params:
            fraction=0.8,
            matches=10,
            prefix="SSU_",
            email=config["ncbi_email"],
            options="Metazoa,Viridiplantae"
        shell:
            r"""
            python {SCRIPTS_DIR}/filter_reads_assign_taxonomy.py \
                --id_info {input.id_info} \
                --fraction {params.fraction} \
                --matches {params.matches} \
                --prefix {params.prefix} \
                --out {output.tsv} \
                --read_id {output.read_ids} \
                --failed_id {output.failed_ids} \
                --discarded_id {output.discarded_ids} \
                --taxonomy_discarded_id {output.taxonomy_discarded_id} \
                --statistics {output.statistics} \
                --discard {params.options} \
                --email {params.email} \
                --taxonomy {input.taxonomy} &> {log.stdout}
            """
        
    rule rrna_extraction:
        conda:
            "../envs/seqkit.yaml"
        message:
            "[Seqkit] extract rRNA reads for {wildcards.sample}"
        input:
            cleaned_r1=f"{RESULTS_DIR}/qc/{{sample}}/decontamination/{{sample}}_R1.cleaned.fastq.gz",
            cleaned_r2=f"{RESULTS_DIR}/qc/{{sample}}/decontamination/{{sample}}_R2.cleaned.fastq.gz",
            read_ids=rules.filter_and_assign_taxonomy.output.read_ids
        output:
            rrna_r1=f"{RESULTS_DIR}/rRNA/{{sample}}/filtered/{{sample}}_rRNA_1.fastq.gz",
            rrna_r2=f"{RESULTS_DIR}/rRNA/{{sample}}/filtered/{{sample}}_rRNA_2.fastq.gz"
        log:
            stdout=f"{RESULTS_DIR}/rRNA/{{sample}}/logs/rrna_extraction.log"
        benchmark:
            f"{RESULTS_DIR}/rRNA/{{sample}}/benchmarks/{{sample}}_rrna_extraction.txt"
        threads:
            config["qc"]["decontamination"]["threads"]
        params:
            prefix=f"{RESULTS_DIR}/rRNA/{{sample}}/filtered/{{sample}}_rRNA"
        shell:
            r"""
            set -euo pipefail

            mkdir -p $(dirname {output.rrna_r1})

            seqkit grep \
                -f {input.read_ids} \
                --paired \
                --threads {threads} \
                -1 {input.cleaned_r1} \
                -2 {input.cleaned_r2} \
                -o {params.prefix} \
                &> {log.stdout}
            """
    
    rule non_rrna_extraction:
        conda:
            "../envs/seqkit.yaml"
        message:
            "[Seqkit] extract non-rRNA reads for {wildcards.sample}"
        input:
            cleaned_r1=f"{RESULTS_DIR}/qc/{{sample}}/decontamination/{{sample}}_R1.cleaned.fastq.gz",
            cleaned_r2=f"{RESULTS_DIR}/qc/{{sample}}/decontamination/{{sample}}_R2.cleaned.fastq.gz",
            readid_txt=rules.extract_all_aligned_ids.output.readid_txt
        output:
            non_rrna_r1=f"{RESULTS_DIR}/nonrRNA/{{sample}}/filtered/{{sample}}_nonrRNA_1.fastq.gz",
            non_rrna_r2=f"{RESULTS_DIR}/nonrRNA/{{sample}}/filtered/{{sample}}_nonrRNA_2.fastq.gz"
        log:
            stdout=f"{RESULTS_DIR}/nonrRNA/{{sample}}/logs/non_rrna_extraction.log"
        benchmark:
            f"{RESULTS_DIR}/nonrRNA/{{sample}}/benchmarks/{{sample}}_non_rrna_extraction.txt"
        threads:
            config["qc"]["decontamination"]["threads"]
        params:
            prefix=f"{RESULTS_DIR}/nonrRNA/{{sample}}/filtered/{{sample}}_nonrRNA"
        shell:
            r"""
            set -euo pipefail

            mkdir -p $(dirname {output.non_rrna_r1})

            seqkit grep \
                -v -f {input.readid_txt} \
                --paired \
                --threads {threads} \
                -1 {input.cleaned_r1} \
                -2 {input.cleaned_r2} \
                -o {params.prefix} \
                &> {log.stdout}
            """

