#!/bin/bash

t=10
minimap2 -t "$t" -a ref_genome.fa sample_corrected.renamed.fasta > sample.sam
samtools view -@ "$t" -bo sample.bam sample.sam
samtools sort -@ "$t" -o sorted_sample.bam sample.bam
samtools index -@ "$t" sorted_sample.bam
