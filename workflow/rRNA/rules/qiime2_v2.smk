#https://docs.qiime2.org/2024.10/tutorials/importing/ -> “Fastq manifest” formats

rule qiime_manifest:
    message:
        "[QIIME2] generate manifest file across samples for stage {stage}"
    input:
        cleaned_r1=f"{RESULTS_DIR}/{{sample}}_R1.cleaned.fastq.gz", #insert new filtered reads
    output:
        manifest="{RESULTS_DIR}/rRNA/qiime/{stage}/manifest.csv"
                aligned=f"{RESULTS_DIR}/RNA/{{sample}}/bbmap/{{sample}}_aligned.sam",

    run:
        import os

        os.makedirs(os.path.dirname(output.manifest), exist_ok=True)

        with open(output.manifest, "w") as fh:
            fh.write(
                "sample-id\tabsolute-filepath\n"
            )
                echo -e "${sample}\t$(realpath "$fwd")" >> "$OUTPUT_MANIFEST"


            for sample, fwd in zip(SAMPLES, input.cleaned_r1):
                fh.write(
                    f"{sample}\t{os.path.abspath(fwd)}\n"
                )

#https://amplicon-docs.qiime2.org/en/stable/how-to-guides/how-to-import/
# import your demultiplexed sequences,
rule qiime_import:
    conda: "../envs/qiime2.yaml"
    message:
        "[QIIME2] import reads for {wildcards.sample} for stage {wildcards.stage}"
    input:
        manifest=rules.qiime_manifest.output.manifest
    output:
        qza="{outdir}/{sample}/rRNA/qiime/{stage}/{sample}_{stage}_demux.qza"
    log:
        stdout = "{outdir}/{sample}/logs/qiime_import.log"
    benchmark:
        "{outdir}/{sample}/benchmarks/qiime_import.txt"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.qza})

        qiime tools import \
          --type 'SampleData[PairedEndSequencesWithQuality]' \
          --input-format PairedEndFastqManifestPhred33V2 \
          --input-path {input.manifest} \
          --output-path {output.qza} \
          > {log.stdout} 2>&1
        """

#qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path manifest.tsv \
  --output-path 01_qiime_import/demux.qza \
  --input-format SingleEndFastqManifestPhred33V2
