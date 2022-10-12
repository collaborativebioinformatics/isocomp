
## File spec for the unique isoforms

### Header

Do we need a header? We should return to this later....

Ideas:

```
##Unique isoforms
##Samples: 
##Reference: 
##CHR    WINDOW_START    WINDOW_END    CHR    SHARED    ISOFORM_START    ISOFORM_END ...
```

The samples should be a comma-separated string of the samples

The reference should include the reference name for alignments

We could also include things like the alignment algorithm + version, etc. 

### Columns

**1. CHR**

The chromosome. This either includes a prefix 'chr' or not. String.

e.g. `chr1`, `3`, `X`, `chrX`


**2. WINDOW_START**

This is the coorindate of where our bin/window begins. Integer.
It is a 1-based index location on the chromosome.

(I'm happy with the name `BIN_START` as well...)


**3. WINDOW_END**

This is the coorindate of where our bin/window begins. Integer.
It is a 1-based index location on the chromosome.

(I'm happy with the name `BIN_END` as well...)



**4. SHARED_PAIRS**

Based on the sample names given in the algorithm in the FASTA files (?), 
this should be a comma-separated string of pairs of samples whereby the unique isoform listed is found

String.

e.g.
We have three samples: s1, s2, s3

There's a unique isoform in s2, i. we will find this unique isoform is unique against s1 and s2 when calculating intersections. 

We could either output the pairs

`s2_s3,s2_s1`

Otherwise, we could write out the sample where the unique isoforms are found

```
s2
```

But this will require an extra step of string matching. 

e.g. consider the case if there's a unique isoform in s1 and s2. It could be the same isoform, or a different isoform. 



**5. ISOFORM_START**

Start of the alignment of the isoform (against T2T). We get this from the BAM.


**6. ISOFORM_END**

End of the alignment of the isoform (against T2T). We get this from the BAM.

**6. ISOFORM_NAME**

QNAME of the isoform from the alignments (which were originally given in the FASTA file)


**7. TXNAME

The transcript name associated with this isoform

**8. GENE_NAME

The common gene name (based on HGNC) for the gene. It would be best if these were gene names, not Ensembl IDs.

**9. SEQUENCE

The entire isoform sequence.












