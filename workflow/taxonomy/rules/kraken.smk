rule kraken2:
    conda:
        "../envs/kraken2.yaml"
    wildcard_constraints:
        marker="SSU|ITS"
    message:
        "[Kraken2] Classifying {wildcards.marker} reads from {wildcards.sample} against {wildcards.database}"
    input:
        read_1=lambda wc: taxonomy_read_1(
            wc,
            results_dir=RESULTS_DIR,
        ),
        read_2=lambda wc: taxonomy_read_2(
            wc,
            results_dir=RESULTS_DIR,
        ),
    output:
        report=f"{TAXONOMY_DIR}/{{sample}}/{{marker}}/{{database}}/kraken/{{sample}}.{{database}}.k2report",
        kraken=f"{TAXONOMY_DIR}/{{sample}}/{{marker}}/{{database}}/kraken/{{sample}}.{{database}}.kraken",
    log:
        stdout=f"{TAXONOMY_DIR}/{{sample}}/{{marker}}/{{database}}/logs/{{sample}}.{{database}}.kraken2.log",
    benchmark:
        f"{TAXONOMY_DIR}/{{sample}}/{{marker}}/{{database}}/benchmarks/{{sample}}.{{database}}.kraken2.txt",
    params:
        db=lambda wc: kraken_db(wc, config=config),
        options=lambda wc: kraken_options(wc, config=config),
    threads:
        lambda wc: kraken_threads(wc, config=config)
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.kraken}) $(dirname {log.stdout})

        r1_lines=$(gzip -cd {input.read_1} | wc -l)
        r2_lines=$(gzip -cd {input.read_2} | wc -l)
        
        if [ "$r1_lines" -eq 0 ] && [ "$r2_lines" -eq 0 ]; then
            echo "[Kraken2] No input read pairs; skipping classification." > {log.stdout}
            : > {output.kraken}
            : > {output.report}
        elif [ "$r1_lines" -ne "$r2_lines" ] || [ "$((r1_lines % 4))" -ne 0 ]; then
            echo "[ERROR] Paired FASTQs have unequal line counts or incomplete records." > {log.stdout}
            exit 1
        else
            kraken2 \
                --db {params.db} \
                --threads {threads} \
                --report {output.report} \
                --output {output.kraken} \
                {input.read_1} \
                {input.read_2} \
                {params.options} \
                > {log.stdout} 2>&1
            if [ ! -s {output.kraken} ] || [ ! -s {output.report} ]; then
                echo "[ERROR] Missing or empty Kraken2 outputs for nonempty inputs." >> {log.stdout}
                exit 1
            fi
        fi
        """

rule kraken2_group_reads:
    conda:
        "../envs/krakentools.yaml"
    message:
        "[KrakenTools] Separating SILVA-classified SSU reads from {wildcards.sample} into prokaryotes and eukaryotes"
    input:
        read_1=lambda wc: original_ssu_read(wc, 1, results_dir=RESULTS_DIR),
        read_2=lambda wc: original_ssu_read(wc, 2, results_dir=RESULTS_DIR),
        report=f"{TAXONOMY_DIR}/{{sample}}/SSU/silva/kraken/{{sample}}.silva.k2report",
        kraken=f"{TAXONOMY_DIR}/{{sample}}/SSU/silva/kraken/{{sample}}.silva.kraken"
    output:
        prokaryotes_1=grouped_ssu_read("{sample}", "prokaryotes", 1, results_dir=RESULTS_DIR),
        prokaryotes_2=grouped_ssu_read("{sample}", "prokaryotes", 2, results_dir=RESULTS_DIR),
        eukaryotes_1=grouped_ssu_read("{sample}", "eukaryotes", 1, results_dir=RESULTS_DIR),
        eukaryotes_2=grouped_ssu_read("{sample}", "eukaryotes", 2, results_dir=RESULTS_DIR)
    log:
        stdout=f"{TAXONOMY_DIR}/{{sample}}/SSU/silva/logs/{{sample}}.group_reads.log"
    params:
        prokaryote_taxids=" ".join(map(str, config["taxonomy"]["SSU"]["domain_split"]["prokaryotes"])),
        eukaryote_taxids=" ".join(map(str, config["taxonomy"]["SSU"]["domain_split"]["eukaryotes"]))
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.prokaryotes_1}) $(dirname {output.eukaryotes_1}) $(dirname {log.stdout})
        (
            if [ ! -s {input.kraken} ] && [ ! -s {input.report} ]; then
                echo "[KrakenTools] Empty SILVA classification; creating empty read groups."
                : > {output.prokaryotes_1}.tmp.fastq
                : > {output.prokaryotes_2}.tmp.fastq
                : > {output.eukaryotes_1}.tmp.fastq
                : > {output.eukaryotes_2}.tmp.fastq
            elif [ ! -s {input.kraken} ] || [ ! -s {input.report} ]; then
                echo "[ERROR] Inconsistent SILVA classification outputs: only one file is nonempty."
                exit 1
            else
                extract_kraken_reads.py \
                    -k {input.kraken} \
                    -r {input.report} \
                    -s1 {input.read_1} \
                    -s2 {input.read_2} \
                    -t {params.prokaryote_taxids} \
                    --include-children \
                    --fastq-output \
                    -o {output.prokaryotes_1}.tmp.fastq \
                    -o2 {output.prokaryotes_2}.tmp.fastq
                extract_kraken_reads.py \
                    -k {input.kraken} \
                    -r {input.report} \
                    -s1 {input.read_1} \
                    -s2 {input.read_2} \
                    -t {params.eukaryote_taxids} \
                    --include-children \
                    --fastq-output \
                    -o {output.eukaryotes_1}.tmp.fastq \
                    -o2 {output.eukaryotes_2}.tmp.fastq
            fi
            gzip -c {output.prokaryotes_1}.tmp.fastq > {output.prokaryotes_1}
            gzip -c {output.prokaryotes_2}.tmp.fastq > {output.prokaryotes_2}
            gzip -c {output.eukaryotes_1}.tmp.fastq > {output.eukaryotes_1}
            gzip -c {output.eukaryotes_2}.tmp.fastq > {output.eukaryotes_2}
            rm -- \
                {output.prokaryotes_1}.tmp.fastq \
                {output.prokaryotes_2}.tmp.fastq \
                {output.eukaryotes_1}.tmp.fastq \
                {output.eukaryotes_2}.tmp.fastq
        ) > {log.stdout} 2>&1
        """