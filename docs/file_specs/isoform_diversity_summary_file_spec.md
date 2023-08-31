
# Motivation

The main point of this project is to summarize the isoform diversity between
input N samples. 

Along with a file reporting the “variants” between pairwise alignments, we also
need to output a summary file for user. 

That is, we need a summary stating if there is “variation” between two isoforms
(during the pairwise alignment) or there is no “variation”. 


# File specification, same isoforms found

### Header

This lists all of the "same" isoforms we find for the user.

The "unfiltered" file lists the isoforms found which have no mismatches after
doing pairwise sequence alignment. 

We should also calculate a "filtered" file whereby we remove the alignments
which have mismatches solely due to the length differing between the REFERENCE
and QUERY (see "False Positive" discussion elsewhere).

```
##Pairs of same isoforms, unfiltered #Samples: #Reference: #CHR    PAIR_ID
#WINDOW_START    WINDOW_END    ISOFORM1_START    ISOFORM1_END    ISOFORM2_START
#ISOFORM2_END
```

For the "filtered" file, this should be:

```
##Pairs of same isoforms, filtered #Samples: #Reference: #CHR    PAIR_ID
#WINDOW_START    WINDOW_END    ISOFORM1_START    ISOFORM1_END    ISOFORM2_START
#ISOFORM2_END
```



### Columns

**1. CHR**

The chromosome. This either includes a prefix 'chr' or not. String.

e.g. `chr1`, `3`, `X`, `chrX`


**2. PAIR_ID**

An underscore-separated pair of the isoform "names", i.e. the QNAMEs of the
isoform from the alignments (which were originally given in the FASTA file)


Based on the sample names given in the algorithm in the FASTA files (?), this
should be a comma-separated string of pairs of samples 

String.

e.g.  We have three samples: s1, s2, s3

There's a unique isoform in s2, i. we will find this unique isoform is unique
against s1 and s2 when calculating intersections. 

The id will be `s1_s2` 



**3. WINDOW_START**

This is the coorindate of where our bin/window begins. Integer.  It is a 1-based
index location on the chromosome.

(I'm happy with the name `BIN_START` as well...)


**4. WINDOW_END**

This is the coorindate of where our bin/window begins. Integer.  It is a 1-based
index location on the chromosome.

(I'm happy with the name `BIN_END` as well...)



**5. ISOFORM1_START**

Start of the alignment of the first isoform (against T2T). We get this from the
BAM.


**6. ISOFORM1_END**

End of the alignment of the first isoform (against T2T). We get this from the
BAM.



**7. ISOFORM2_START**

Start of the alignment of the second isoform (against T2T). We get this from the
BAM.


**8. ISOFORM2_END**

End of the alignment of the second isoform (against T2T). We get this from the
BAM.

# File specification, different isoforms found

I think that the "variant" file we are writing covers this information....

Do we need a summary file here? Just to summarize for the users that which pairs
of isoforms are different and maybe calculate the mismatch score?

Something like the following:



```
##Pairs of different isoforms, unfiltered #Samples: #Reference: #CHR    PAIR_ID
#WINDOW_START    WINDOW_END    PERCENT_MISMATCH   ISOFORM1_START    ISOFORM1_END
#ISOFORM2_START    ISOFORM2_END    
```


### Columns

Same as above, but there is a `PERCENT_MISMATCH` field.


**PERCENT_MISMATCH**

The percentage mismatch between the two 



