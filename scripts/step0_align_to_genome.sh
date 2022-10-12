#!/bin/bash

"""
For some background information:
  Step 0)

  With three input FASTAs (e.g. HG002, HG004, HG005), for each fasta file

  (0B) Make sure FASTA seq IDs are both unique and have the
  sample name. Proposal: "SampleName_i_existingFastaSeqID"
  whereby i is the increment

  (0C) Align these three FASTAs to the T2T reference, and
  create a BAM file. This script does this part of step 0.
"""

t=10
minimap2 -t "$t" -a ref_genome.fa sample_corrected.renamed.fasta > sample.sam
samtools view -@ "$t" -bo sample.bam sample.sam
samtools sort -@ "$t" -o sorted_sample.bam sample.bam
samtools index -@ "$t" sorted_sample.bam
