import sys
import re
def newnames(oldname, n, kmers=None, overlap=None, header=None):
    """
    >>> newnames('some.fasta', 1)
    ['some.split.fasta']
    >>> newnames('some.fasta', 2)
    ['some.0.fasta', 'some.1.fasta']
    >>> newnames('some', 2)
    ['some.0', 'some.1']
    >>> newnames('some.fasta', 2, kmers=1000)
    ['some.0.1Kmer.fasta', 'some.1.1Kmer.fasta']
    >>> newnames('some.fasta', 2, kmers=10000, overlap=2000)
    ['some.0.10Kmer.2Koverlap.fasta', 'some.1.10Kmer.2Koverlap.fasta']
    >>> newnames('some.fasta', 1, kmers=100000, overlap=2000)
    ['some.split.100Kmer.2Koverlap.fasta']
    """
    if kmers and kmers % 1000 == 0: kmers = "%iK" % (kmers/1000)
    if overlap and overlap % 1000 == 0: overlap = "%iK" % (overlap/1000)

    p = oldname.rfind("fa")
    kstr = kmers is not None and ("%smer." % kmers) or ""
    ostr = overlap is not None and ("%soverlap." % overlap) or ""
    if p != -1:
        pattern = oldname[:p] + "%s." + kstr + ostr + oldname[p:]
    else:
        pattern = oldname + kstr + ostr + ".%s"
    if n == 1:
        names = [pattern % "split"]
    else:
        width = len(str(n))
        names = [pattern % str(i).rjust(width, '0') for i in range(n)]
    return names


rule translate_filtered_trinity:
    input:
        ancient(f"results/mRNA/trinity/contigs_ncrna_filtered.fasta"),
    output:
        f"results/mRNA/trinity/translated/contigs_ncrna_filtered.fasta",
    params:
        orfs=int(config.get("ORFs_translate", 6)),
    conda:
        "../envs/align_contigs_to_database.yaml"
    log:
        "logs/align_contigs_to_database/translate.log",
    benchmark:
        "benchmarks/align_contigs_to_database/translate.log"
    shell:
        "transeq -sformat pearson -clean -frame "
        "{params.orfs} -sequence {input} -outseq {output} "
        ">> {log} 2>&1"


rule split_fasta:
    input:
        f"results/mRNA/trinity/translated/contigs_ncrna_filtered.fasta",
    output:
        newnames(f"results/mRNA/trinity/translated/contigs_ncrna_filtered.fasta", int(config.get("split_fasta", 2)))
    params:
        n=int(config.get("split_fasta", 2)),
    conda:
        "../envs/align_contigs_to_database.yaml"
    log:
        "logs/align_contigs_to_database/split.log",
    benchmark:
        "benchmarks/align_contigs_to_database/split.log"
    shell:
        "pyfasta split -n {params.n} {input} "
        ">> {log} 2>&1"


rule sword:
    input:
        f"results/mRNA/trinity/translated/contigs_ncrna_filtered.{{index}}.fasta",
    output:
        f"results/mRNA/trinity/sword/temp/{{database}}.{{index}}.tsv",
    conda:
        "../envs/align_contigs_to_database.yaml"
    params:
        database_path=lambda wildcards: sword_databases.get(wildcards["database"]),
    log:
        f"logs/align_contigs_to_database/sword_{{database}}_{{index}}.log",
    benchmark:
        f"benchmarks/align_contigs_to_database/sword_{{database}}_{{index}}.log"
    threads: int(config.get("align_contigs_to_database-THREADS", 50))
    shell:

        "sword -i {input} "
        "-t {threads} -o {output} "
        "-f bm9 -j {params.database_path} -c 30000 "
        ">> {log} 2>&1"


rule align_contigs_to_database:
    input:
        expand(
            "results/mRNA/trinity/sword/temp/{{database}}.{index}.tsv",
            index = [re.search('contigs_ncrna_filtered.(.*).fasta', file).group(1) for file in newnames(f"results/mRNA/trinity/translated/contigs_ncrna_filtered.fasta", int(config.get("split_fasta", 2)))]
        ),
    output:
        f"results/mRNA/ML_SWORD_{{database}}_result.tsv",
    conda:
        "../envs/align_contigs_to_database.yaml"
    log:
        f"logs/align_contigs_to_database/cat_{{database}}.log",
    benchmark:
        "benchmarks/align_contigs_to_database/cat_{{database}}.log"
    shell:
        "cat {input} > {output} && "
        "echo Finished > {log}"
