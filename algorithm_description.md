
## Algorithm idea

### Step (0) 
Take three input FASTAs (e.g. HG002, HG004, HG005)
(0B) Make sure FASTA seq IDs are both unique and have the sample name. Proposal: "SampleName_i_existingFastaSeqID" whereby i is the increment
(0C) Align these three FASTAs to the T2T reference, and create a BAM file

### Step (1) 
Calculate coverage in BAM file using a rolling window of N bp whereby N=100.
The output of this coverage file should be a tab-delimited BED-like file for each region whereby the coverage is >0
i.e.
```
START    END   COVERAGE
```

### Step (2) 
Download a "gene regions" BED file with [tss, end] coordinates for each gene


### Step (3) 
Calculate Bedtools intersect between these two BED files. Output the intersect BED files. These windows should have a unique name.

### Step (4)
For each [tss, end] window, calculate the bins starts/ends  (as @Yutong Qiu describes)

**NOTE:** we need to output this file for users. (Why? Because (A) it's useful for computational biologists to see where the original windows where. (B) we don't want to give summaries for windows with no unique isoforms; by default, all isoforms between these samples are shared and the sequences are the exact same.)

The output file should be a simple BED file:

````
window_name window_start window_end
````

### Step (5)
For each window (i.e. using the start and end coordinate for each window), extract the "reads" i.e. the isoforms which aligned in that window. 

How? We first must extract the QNAMEs for any isoform sequence which aligned in that window. 

Based on these QNAMEs, we can categorize these into the different samples, i.e. QNAMEs for sample1, QNAMEs for sample2, QNAMEs for sample3. 

We extract the full isoform sequences for each of these QNAMEs. 

One idea:
We can output this into a separate FASTA file for each window. We'll need to name the output FASTA with the unique name given in Step (3)

### Step (6)
We now have sequences from three samples (HG002, HG004, HG005) for each window.

We will try doing intersections of sets first, based on the number of combinations.
e.g. "sample1 vs sample2", "sample2 vs sample 3", "sample1 vs sample3"

i.e. 
```
S1 + S2
S2 + S3
S1 + S3
```

If the intersection between these sequences equals the two samples, there is nothing "unique" between these pairs. 

If there is something special, we have a unique isoform in that window. 






