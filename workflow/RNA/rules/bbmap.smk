rule bbmap:
    conda:
        "../envs/bbmap.yaml"
    message:
        "[BBMap] Running BBMap for {wildcards.sample}"
    input:
        ribodetector_r1 = f"{RESULTS_DIR}/RNA/{{sample}}/ribodetector/{{sample}}.rRNA.r1.fastq.gz",
        ribodetector_r2 = f"{RESULTS_DIR}/RNA/{{sample}}/ribodetector/{{sample}}.rRNA.r2.fastq.gz",
        ref=config["databases"]["sortmeRNA_ssu"]
    output:
        aligned=f"{RESULTS_DIR}/RNA/{{sample}}/bbmap/{{sample}}_aligned.bam",
        unaligned=f"{RESULTS_DIR}/RNA/{{sample}}/bbmap/{{sample}}_unaligned.bam",
    log:
        stdout = f"{RESULTS_DIR}/RNA/{{sample}}/logs/bbmap.log"
    benchmark:
        f"{RESULTS_DIR}/RNA/{{sample}}/benchmarks/bbmap.txt"
    params:
        options=config["RNA"]["bbmap"]["options"]
    threads:
        config["RNA"]["bbmap"]["threads"]
    shell:
        r"""
        set -euo pipefail
        
        bbmap.sh \
            in1={input.ribodetector_r1} \
            in2={input.ribodetector_r2} \
            ref={input.ref} \
            outm={output.aligned} \
            outu={output.unaligned} \
            threads={threads} \
            {params.options} \
            > {log.stdout} 2>&1
        """

rule rRNA_id:
    conda:
        "../envs/bbmap.yaml"
    message:
        "[BBMap] identify aligned rRNA reads for {wildcards.sample}"
    input:
        aligned=rules.bbmap.output.aligned,
    output:
        rna_seqid=f"{RESULTS_DIR}/RNA/{{sample}}/bbmap/{{sample}}_rRNA_readid.txt",
    log:
        stdout=f"{RESULTS_DIR}/RNA/{{sample}}/logs/bbmap_seqread_ids.log"
    benchmark:
        f"{RESULTS_DIR}/RNA/{{sample}}/benchmarks/bbmap_seqread_ids.log"
    params:
        filters="-F 4 -q 30" # consider adding additional filters to ensure we get high confidence rRNA reads, such as higher mapping quality, or specific flags to ensure we only get primary alignments, and not secondary or supplementary alignments.
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.rna_seqid})

        samtools view {params.filters} {input.aligned} | cut -f1 | sort | uniq > {output.rna_seqid} 2> {log.stdout}
        """

rule filter_reads_by_id:
    conda:
        "../envs/bbmap.yaml"
    message:
        "[BBMap] filter aligned rRNA reads for {wildcards.sample}"
    input:
        ribodetector_r1 = f"{RESULTS_DIR}/RNA/{{sample}}/ribodetector/{{sample}}.rRNA.r1.fastq.gz",
        ribodetector_r2 = f"{RESULTS_DIR}/RNA/{{sample}}/ribodetector/{{sample}}.rRNA.r2.fastq.gz",
        read_ids=rules.rRNA_id.output.rna_seqid
    output:
        filtered_r1=f"{RESULTS_DIR}/rRNA/{{sample}}/filtered/{{sample}}_R1.rRNA.fastq.gz",
        filtered_r2=f"{RESULTS_DIR}/rRNA/{{sample}}/filtered/{{sample}}_R2.rRNA.fastq.gz",
    log:
        stdout=f"{RESULTS_DIR}/rRNA/{{sample}}/logs/filtered.log"
    benchmark:
        f"{RESULTS_DIR}/rRNA/{{sample}}/benchmarks/{{sample}}_filtered.txt"
    threads:
        config["RNA"]["bbmap"]["threads"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.filtered_r1})

        seqkit grep \
            -f {input.read_ids} \
            -v \
            --threads {threads} \
            {input.ribodetector_r1} \
            -o {output.filtered_r1} \
            2> {log.stdout}

        seqkit grep \
            -f {input.read_ids} \
            -v \
            --threads {threads} \
            {input.ribodetector_r2} \
            -o {output.filtered_r2} \
            2>> {log.stdout}
        """