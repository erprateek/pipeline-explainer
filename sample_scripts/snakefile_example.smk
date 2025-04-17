rule all:
    input:
        "results/aligned.bam"

rule fastqc:
    input:
        "data/sample.fastq"
    output:
        "qc/sample_fastqc.html"
    shell:
        "fastqc {input} -o qc/"

rule align_reads:
    input:
        "data/sample.fastq"
    output:
        "results/aligned.bam"
    shell:
        "bwa mem ref.fa {input} | samtools view -Sb - > {output}"

