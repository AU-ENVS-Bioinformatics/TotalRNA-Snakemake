import gzip
import statistics

def median_read_length(fq, n_reads=10000):
    lengths = []

    with gzip.open(fq, "rt") as handle:
        for i, line in enumerate(handle):
            # sequence line in FASTQ
            if i % 4 == 1:
                lengths.append(len(line.strip()))

                if len(lengths) >= n_reads:
                    break

    return round(statistics.median(lengths)) if lengths else 150

def nearest_read_model(length):
    allowed = [50, 75, 100, 125, 150, 200, 250]
    return min(allowed, key=lambda x: abs(x - length))

def ribodetector_len(wc):
    fq = (
        RESULTS_DIR
        / "qc"
        / wc.sample
        / "decontamination"
        / f"{wc.sample}_R1.cleaned.fastq.gz"
    )

    median_len = median_read_length(fq)

    model_len = nearest_read_model(median_len)

    print(
        f"[RiboDetector] {wc.sample}: "
        f"median read length = {median_len}, "
        f"using model length = {model_len}"
    )

    return model_len

rule ribodetector:
    conda:
        "../envs/ribodetector.yaml"
    message:
        "[RiboDetector] Running RiboDetector for {wildcards.sample}"
    input:
        cleaned_r1=f"{RESULTS_DIR}/qc/{{sample}}/decontamination/{{sample}}_R1.cleaned.fastq.gz",
        cleaned_r2=f"{RESULTS_DIR}/qc/{{sample}}/decontamination/{{sample}}_R2.cleaned.fastq.gz",
    output:
        nonrna_r1=f"{RESULTS_DIR}/RNA/{{sample}}/ribodetector/{{sample}}_nonrRNA_1.fastq.gz",
        nonrna_r2=f"{RESULTS_DIR}/RNA/{{sample}}/ribodetector/{{sample}}_nonrRNA_2.fastq.gz",
        rna_r1=f"{RESULTS_DIR}/RNA/{{sample}}/ribodetector/{{sample}}_rRNA_1.fastq.gz",
        rna_r2=f"{RESULTS_DIR}/RNA/{{sample}}/ribodetector/{{sample}}_rRNA_2.fastq.gz",
    log:
        stdout = f"{RESULTS_DIR}/RNA/{{sample}}/logs/ribodetector.log"
    benchmark:
        f"{RESULTS_DIR}/RNA/{{sample}}/benchmarks/ribodetector.txt"
    params:
        options=config["RNA"]["ribodetector"]["options"],
        read_len=ribodetector_len
    threads:
        config["RNA"]["ribodetector"]["threads"]
    shell:
        r"""
        set -euo pipefail
        
        ribodetector_cpu \
            --input {input.cleaned_r1} {input.cleaned_r2} \
            --output {output.nonrna_r1} {output.nonrna_r2} \
            --rrna {output.rna_r1} {output.rna_r2} \
            --threads {threads} \
            --len {params.read_len} \
            {params.options} \
            --log {log.stdout} 2>&1
        """

rule link_non_rRNA_ribodetector:
    message:
        "[RiboDetector] Linking non-rRNA reads for {wildcards.sample}"
    input:
        nonrna_r1=rules.ribodetector.output.nonrna_r1,
        nonrna_r2=rules.ribodetector.output.nonrna_r2,
    output:
        linked_r1=f"{RESULTS_DIR}/nonrRNA/{{sample}}/filtered/{{sample}}_nonrRNA_1.fastq.gz",
        linked_r2=f"{RESULTS_DIR}/nonrRNA/{{sample}}/filtered/{{sample}}_nonrRNA_2.fastq.gz",
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.linked_r1})

        ln -s {input.nonrna_r1} {output.linked_r1}
        ln -s {input.nonrna_r2} {output.linked_r2}
        """