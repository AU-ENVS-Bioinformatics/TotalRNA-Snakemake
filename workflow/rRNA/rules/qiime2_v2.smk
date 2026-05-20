#https://docs.qiime2.org/2024.10/tutorials/importing/ -> “Fastq manifest” formats

rule qiime_manifest:
    message:
        "[QIIME2] generate manifest file across samples for stage {stage}"

    input:
        fwd=lambda wc: [
            f"{SAMPLE_OUTDIR[s]}/{s}/QC/sortmerna/{wc.stage}/{s}_{wc.stage}.aligned_fwd.fq.gz"
            for s in SAMPLES
        ],
        rev=lambda wc: [
            f"{SAMPLE_OUTDIR[s]}/{s}/QC/sortmerna/{wc.stage}/{s}_{wc.stage}.aligned_rev.fq.gz"
            for s in SAMPLES
        ]

    output:
        manifest="{outdir}/rRNA/qiime/{stage}/manifest.csv"

    run:
        import os

        os.makedirs(os.path.dirname(output.manifest), exist_ok=True)

        with open(output.manifest, "w") as fh:
            fh.write(
                "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n"
            )

            for sample, fwd, rev in zip(SAMPLES, input.fwd, input.rev):
                fh.write(
                    f"{sample}\t{os.path.abspath(fwd)}\t{os.path.abspath(rev)}\n"
                )

#https://amplicon-docs.qiime2.org/en/stable/how-to-guides/how-to-import/
# import your demultiplexed sequences,
rule qiime_import:
    conda: "../envs/qiime2.yaml"
    message:
        "[QIIME2] import reads for {wildcards.sample} for stage {wildcards.stage}"
    input:
        manifest="{outdir}/rRNA/qiime/{stage}/manifest.csv"
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

#single merged amplicon sequences
rule qiime_join_pairs:
    conda: "../envs/qiime2.yaml"
    message:
        "[QIIME2] join paired-end reads for {wildcards.sample} for stage {wildcards.stage}"
    input:
        demux="{outdir}/{sample}/rRNA/qiime/{stage}/{sample}_{stage}_demux.qza"
    output:
        joined="{outdir}/{sample}/rRNA/qiime/{stage}/{sample}_{stage}_joined.qza"
    log:
        stdout = "{outdir}/{sample}/logs/qiime_join_pairs.log"
    benchmark:
        "{outdir}/{sample}/benchmarks/qiime_join_pairs.txt"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.joined})

        qiime vsearch join-pairs \
          --i-demultiplexed-seqs {input.demux} \
          --o-joined-sequences {output.joined}
        """

rule qiime_quality_filter:
    conda: "../envs/qiime2.yaml"
    message:
        "[QIIME2] filter quality scores for {wildcards.sample} for stage {wildcards.stage}"
    input:
        joined="{outdir}/{sample}/rRNA/qiime/{stage}/{sample}_{stage}_joined.qza"
    output:
        filtered="{outdir}/{sample}/rRNA/qiime/{stage}/{sample}_{stage}_filtered.qza",
        stats="{outdir}/{sample}/rRNA/qiime/{stage}/{sample}_{stage}_filter-stats.qza"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.filtered})

        qiime quality-filter q-score \
          --i-demux {input.joined} \
          --o-filtered-sequences {output.filtered} \
          --o-filter-stats {output.stats}
        """

# Begin OTU creation steps
rule qiime_dereplicate:
    conda: "../envs/qiime2.yaml"
    message:
        "[QIIME2] dereplicate sequences for {wildcards.sample} for stage {wildcards.stage}"
    input:
        filtered="{outdir}/{sample}/rRNA/qiime/{stage}/{sample}_{stage}_filtered.qza"
    output:
        table="{outdir}/{sample}/rRNA/qiime/{stage}/{sample}_{stage}_table-derep.qza",
        seqs="{outdir}/{sample}/rRNA/qiime/{stage}/{sample}_{stage}_seqs-derep.qza"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.table})

        qiime vsearch dereplicate-sequences \
          --i-sequences {input.filtered} \
          --o-dereplicated-table {output.table} \
          --o-dereplicated-sequences {output.seqs}
        """

#######################################
#### QIIME2 de-novo OTU clustering ####
#######################################
#   with vsearch at 97% identity, which is the default for OTU clustering.
#   OTUs defined only by within‑dataset similarity
#   Ideal for:
#       novel environments
#       non‑standard genes
#       mixed or incomplete references

rule qiime_cluster_denovo:
    conda: "../envs/qiime2.yaml"
    message:
        "[QIIME2] cluster de-novo OTUs for {wildcards.sample} for stage {wildcards.stage}"
    input:
        table="{outdir}/{sample}/rRNA/qiime/{stage}/{sample}_{stage}_table-derep.qza",
        seqs="{outdir}/{sample}/rRNA/qiime/{stage}/{sample}_{stage}_seqs-derep.qza"
    output:
        otu_table="{outdir}/{sample}/rRNA/qiime/{stage}/{sample}_{stage}_otu-table.qza",
        otu_seqs="{outdir}/{sample}/rRNA/qiime/{stage}/{sample}_{stage}_otu-seqs.qza"
    params:
        id=0.97
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.otu_table})

        qiime vsearch cluster-features-de-novo \
          --i-table {input.table} \
          --i-sequences {input.seqs} \
          --p-perc-identity {params.id} \
          --o-clustered-table {output.otu_table} \
          --o-clustered-sequences {output.otu_seqs}
        """

##################################################
#### Closed‑reference OTUs (reference‑locked) ####
##################################################

#   OTUs must match reference sequences
#   No new OTUs are created
#   Reads without hits are discarded

# But here each OTU is an reference sequence so this makes reproducibility across studies easier, and also allows for taxonomic assignment of the OTUs, but it is not ideal for novel environments or non-standard genes. NB! remember the OTU might still change depending on the sequence data

rule qiime_cluster_closed_ref:
    conda: "../envs/qiime2.yaml"
    message:
        "[QIIME2] cluster closed-reference OTUs for {wildcards.sample} for stage {wildcards.stage}"
    input:
        table="{outdir}/{sample}/rRNA/qiime/{stage}/{sample}_{stage}_table-derep.qza",
        seqs="{outdir}/{sample}/rRNA/qiime/{stage}/{sample}_{stage}_seqs-derep.qza",
        ref="db/{stage}/ref-seqs.qza"
    output:
        otu_table="{outdir}/{sample}/rRNA/qiime/{stage}/{sample}_{stage}_otu-table.qza",
        otu_seqs="{outdir}/{sample}/rRNA/qiime/{stage}/{sample}_{stage}_otu-seqs.qza",
        unmatched="{outdir}/{sample}/rRNA/qiime/{stage}/{sample}_{stage}_unmatched.qza"
    params:
        id=0.97
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.otu_table})

        qiime vsearch cluster-features-closed-reference \
          --i-table {input.table} \
          --i-sequences {input.seqs} \
          --i-reference-sequences {input.ref} \
          --p-perc-identity {params.id} \
          --o-clustered-table {output.otu_table} \
          --o-clustered-sequences {output.otu_seqs} \
          --o-unmatched-sequences {output.unmatched}
        """

######################################
#### Open-reference OTUs (hybrid) ####
######################################
#   OTUs are clustered against reference sequences, but new OTUs are created for reads that do not match the reference
#   Ideal for:
#       novel environments
#       non‑standard genes
#   Provides a balance between the advantages and disadvantages of de-novo and closed-reference clustering.
#   OTUs that match the reference are more likely to be reproducible across studies, while new OTUs can still be created for novel sequences.

rule qiime_cluster_open_ref:
    conda: "../envs/qiime2.yaml"
    message:
        "[QIIME2] cluster open-reference OTUs for {wildcards.sample} for stage {wildcards.stage}"
    input:
        table="{outdir}/{sample}/rRNA/qiime/{stage}/{sample}_{stage}_table-derep.qza",
        seqs="{outdir}/{sample}/rRNA/qiime/{stage}/{sample}_{stage}_seqs-derep.qza",
        ref="db/{stage}/ref-seqs.qza"
    output:
        otu_table="{outdir}/{sample}/rRNA/qiime/{stage}/{sample}_{stage}_otu-table.qza",
        otu_seqs="{outdir}/{sample}/rRNA/qiime/{stage}/{sample}_{stage}_otu-seqs.qza"
    params:
        id=0.97
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.otu_table})

        qiime vsearch cluster-features-open-reference \
          --i-table {input.table} \
          --i-sequences {input.seqs} \
          --i-reference-sequences {input.ref} \
          --p-perc-identity {params.id} \
          --o-clustered-table {output.otu_table} \
          --o-clustered-sequences {output.otu_seqs}
        """

# Chimera filtering of OTUs, which can be done after any of the clustering methods, but is most commonly done after closed-reference clustering, as it is more likely to contain chimeras due to the fact that it is more likely to contain novel sequences that do not match the reference.

rule qiime_chimera_filter:
    conda: "../envs/qiime2.yaml"
    message:
        "[QIIME2] filter chimeras from OTUs for {wildcards.sample} for stage {wildcards.stage}"
    input:
        seqs="{outdir}/{sample}/rRNA/qiime/{stage}/{sample}_{stage}_otu-seqs.qza",
        table="{outdir}/{sample}/rRNA/qiime/{stage}/{sample}_{stage}_otu-table.qza"
    output:
        clean_seqs="{outdir}/{sample}/rRNA/qiime/{stage}/{sample}_{stage}_otu-seqs-nochim.qza",
        chimera_seqs="{outdir}/{sample}/rRNA/qiime/{stage}/{sample}_{stage}_otu-seqs-chim.qza",
        stats="{outdir}/{sample}/rRNA/qiime/{stage}/{sample}_{stage}_chimera-stats.qza",
        clean_table="{outdir}/{sample}/rRNA/qiime/{stage}/{sample}_{stage}_otu-table-nochim.qza"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.clean_table})

        qiime vsearch uchime-denovo \
          --i-sequences {input.seqs} \
          --o-nonchimeras {output.clean_seqs} \
          --o-chimeras {output.chimera_seqs} \
          --o-stats {output.stats}

        qiime feature-table filter-features \
          --i-table {input.table} \
          --m-metadata-file {output.clean_seqs} \
          --o-filtered-table {output.clean_table}
        """

rule qiime_taxonomy:
    conda: "../envs/qiime2.yaml"
    message:
        "[QIIME2] assign taxonomy for {wildcards.sample} for stage {wildcards.stage}"
    input:
        otu_seqs="{outdir}/{sample}/rRNA/qiime/{stage}/{sample}_{stage}_otu-seqs-nochim.qza",
        classifier="db/{stage}/silva_classifier.qza"
    output:
        taxonomy="{outdir}/{sample}/rRNA/qiime/{stage}/{sample}_{stage}_taxonomy.qza"
    log:
        stdout="{outdir}/{sample}/logs/qiime_taxonomy.log"
    benchmark:
        "{outdir}/{sample}/benchmarks/qiime_taxonomy.txt"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.taxonomy})

        qiime feature-classifier classify-sklearn \
          --i-classifier {input.classifier} \
          --i-reads {input.otu_seqs} \
          --o-classification {output.taxonomy} \
          > {log.stdout} 2>&1
        """

