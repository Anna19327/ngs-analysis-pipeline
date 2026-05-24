# =============================================================================
# NGS Alignment Pipeline (Snakemake)
# Tools: fastp → bwa mem → samtools sort/index
# =============================================================================

import os

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

SAMPLES, = glob_wildcards("data/{sample}_R1.fastq.gz")

REF = "reference/genome.fa"

# ---------------------------------------------------------------------------
# Targets
# ---------------------------------------------------------------------------

rule all:
    input:
        expand("results/mapped/{sample}.sorted.bam.bai", sample=SAMPLES),
        expand("results/qc/{sample}_fastp.json", sample=SAMPLES),


# ---------------------------------------------------------------------------
# QC & trimming
# ---------------------------------------------------------------------------

rule fastp:
    input:
        r1 = "data/{sample}_R1.fastq.gz",
        r2 = "data/{sample}_R2.fastq.gz",
    output:
        r1   = "results/trimmed/{sample}_trimmed_R1.fastq.gz",
        r2   = "results/trimmed/{sample}_trimmed_R2.fastq.gz",
        json = "results/qc/{sample}_fastp.json",
        html = "results/qc/{sample}_fastp.html",
    log:
        "logs/fastp/{sample}.log"
    threads: 4
    shell:
        """
        fastp \
            -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            -j {output.json} -h {output.html} \
            -w {threads} \
            2> {log}
        """

# ---------------------------------------------------------------------------
# Alignment
# ---------------------------------------------------------------------------

rule bwa_mem:
    input:
        ref = REF,
        r1  = "results/trimmed/{sample}_trimmed_R1.fastq.gz",
        r2  = "results/trimmed/{sample}_trimmed_R2.fastq.gz",
    output:
        bam = temp("results/mapped/{sample}.bam"),
    params:
        rg = r"@RG\tID:{sample}\tSM:{sample}\tLB:{sample}\tPL:ILLUMINA",
    log:
        "logs/bwa_mem/{sample}.log"
    threads: 8
    shell:
        """
        bwa mem \
            -t {threads} \
            -R '{params.rg}' \
            {input.ref} {input.r1} {input.r2} \
            2> {log} | samtools view -Sb - > {output.bam}
        """

# ---------------------------------------------------------------------------
# Sorting & indexing
# ---------------------------------------------------------------------------

rule samtools_sort:
    input:
        "results/mapped/{sample}.bam",
    output:
        "results/mapped/{sample}.sorted.bam",
    log:
        "logs/samtools_sort/{sample}.log"
    threads: 4
    shell:
        "samtools sort -@ {threads} {input} -o {output} 2> {log}"


rule samtools_index:
    input:
        "results/mapped/{sample}.sorted.bam",
    output:
        "results/mapped/{sample}.sorted.bam.bai",
    log:
        "logs/samtools_index/{sample}.log"
    shell:
        "samtools index {input} 2> {log}"
