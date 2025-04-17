#!/bin/bash

# Step 1: Quality check
fastqc input.fastq -o qc_results/

# Step 2: Alignment
bwa mem -t 8 ref.fa input.fastq > aligned.sam

# Step 3: Convert SAM to BAM
samtools view -Sb aligned.sam > aligned.bam

