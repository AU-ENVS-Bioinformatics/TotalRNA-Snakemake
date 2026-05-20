rule samtools_aligned_rRNA:
    conda:
        "../envs/samtools.yaml"
    message:
        "[Samtools] filter aligned rRNA reads {wildcards.sample}"
    input:
        rRNA_aligned=f"{RESULTS_DIR}/RNA/{{sample}}/ribodetector/{{sample}}.rRNA.aligned.bam",
        SSU_contig_id=config["databases"]["SILVA_138.1_N99"]["reference_contig_id"]
    output:
        rRNA_aligned_read_id=f"{RESULTS_DIR}/qc/{{sample}}/decontamination/{{sample}}_aligned_contaminants_R1.fastq.gz",
    log:
        stdout=f"{RESULTS_DIR}/qc/{{sample}}/logs/aligned_contaminants.log"
    benchmark:
        f"{RESULTS_DIR}/qc/{{sample}}/benchmarks/{{sample}}_bwa_aligned_contaminants.txt"
    params:
        options=
    threads:
        config["qc"]["bwa_contaminants"]["threads"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.rRNA_aligned_fq})

        samtools view -f 4 {input.rRNA_aligned} | awk -v ref_id={input.SSU_contig_id} '$3 == ref_id {print $1}' - | sort -u > {output.rRNA_aligned_read_id}.txt
        """

#bbmap.sh ref=combined.fa in=reads.fq out=mapped.sam ambig=best
#samtools view -f 4 mapped.bam | awk 'NR==FNR {ssu[$1]=1; next} $3 in ssu {print $1}' SSU_contigs.txt - | sort -u > SSU_read_ids.txt
#seqtk subseq reads.fq SSU_read_ids.txt > SSU_reads.fq