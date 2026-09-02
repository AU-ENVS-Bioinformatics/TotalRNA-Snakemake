rule ribodetector:
    conda:
        "../envs/ribodetector.yaml"
    message:
        "[RiboDetector] Running RiboDetector for {wildcards.sample}"
    input:
        cleaned_r1=f"{QC_DIR}/{{sample}}/decontamination/{{sample}}_R1.cleaned.fastq.gz",
        cleaned_r2=f"{QC_DIR}/{{sample}}/decontamination/{{sample}}_R2.cleaned.fastq.gz",
    output:
        nonrna_r1=f"{RNA_CLASSIFIED_DIR}/{{sample}}/non_rRNA/{{sample}}_non_rRNA_1.fastq.gz",
        nonrna_r2=f"{RNA_CLASSIFIED_DIR}/{{sample}}/non_rRNA/{{sample}}_non_rRNA_2.fastq.gz",
        rna_r1=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/ribodetector/rRNA/{{sample}}_rRNA_1.fastq.gz",
        rna_r2=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/ribodetector/rRNA/{{sample}}_rRNA_2.fastq.gz",
    log:
        stdout=f"{RNA_CLASSIFIED_DIR}/{{sample}}/logs/ribodetector_nonrRNA.log"
    benchmark:
        f"{RNA_CLASSIFIED_DIR}/{{sample}}/benchmarks/ribodetector_nonrRNA.txt"
    params:
        options=config["RNA"]["ribodetector"]["options"],
        read_len=150
    threads:
        config["RNA"]["ribodetector"]["threads"]
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.rna_r1})
        mkdir -p $(dirname {log.stdout})

        ribodetector_cpu \
            --input {input.cleaned_r1} {input.cleaned_r2} \
            --output {output.nonrna_r1} {output.nonrna_r2} \
            --rrna {output.rna_r1} {output.rna_r2} \
            --threads {threads} \
            --len {params.read_len} \
            {params.options} \
            --log {log.stdout} 2>&1
        """

# previously we stored all classified sequences within the intermediate and then linked the relevant into the classified, as seen below.
# but in september 2026 i chose to simply directly store the correctly classified reads into the correct folder to remove need for linking 

#rule link_non_rRNA_ribodetector:
#    message:
#        "[RiboDetector] Linking non-rRNA reads for {wildcards.sample}"
#    input:
#        nonrna_r1=rules.ribodetector.output.nonrna_r1,
#        nonrna_r2=rules.ribodetector.output.nonrna_r2,
#    output:
#        linked_r1=f"{RNA_CLASSIFIED_DIR}/{{sample}}/non_rRNA/{{sample}}_non_rRNA_1.fastq.gz",
#        linked_r2=f"{RNA_CLASSIFIED_DIR}/{{sample}}/non_rRNA/{{sample}}_non_rRNA_2.fastq.gz",
#    shell:
#        r"""
#        set -euo pipefail
#
#        mkdir -p $(dirname {output.linked_r1})
#
#        ln -sfn {input.nonrna_r1} {output.linked_r1}
#        ln -sfn {input.nonrna_r2} {output.linked_r2}
#        """

#rule link_rRNA_ribodetector:
#    message:
#        "[RiboDetector] Linking rRNA reads for {wildcards.sample}"
#    input:
#        rna_r1=rules.ribodetector.output.rna_r1,
#        rna_r2=rules.ribodetector.output.rna_r2,
#    output:
#        linked_r1=f"{RNA_CLASSIFIED_DIR}/{{sample}}/rRNA/{{sample}}_rRNA_1.fastq.gz",
#        linked_r2=f"{RNA_CLASSIFIED_DIR}/{{sample}}/rRNA/{{sample}}_rRNA_2.fastq.gz",
#    shell:
#        r"""
#        set -euo pipefail
#
#        mkdir -p $(dirname {output.linked_r1})
#
#        ln -sfn {input.rna_r1} {output.linked_r1}
#        ln -sfn {input.rna_r2} {output.linked_r2}
#        """