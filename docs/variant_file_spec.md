## File spec for the "variants" between the pairwise alignments of isoforms

### Header

Do we need a header? We should return to this later....

Ideas:

```
##Unique isoforms
##Samples: 
##Reference: 
##CHR    REF_START    POS    ID    REF    ALT    REFERENCE    QUERY 
```

The samples should be a comma-separated string of the samples

The reference should include the reference name for alignments

We could also include things like the alignment algorithm + version, etc. 


**Note:**

We are trying to copy the VCF roughly, as this is what users will expect:

```
CHROM   POS   ID   REF   ALT 
```
https://en.wikipedia.org/wiki/Variant_Call_Format


### Columns

**1. CHR**

The chromosome. This either includes a prefix 'chr' or not. String.

e.g. `chr1`, `3`, `X`, `chrX`


**2. REF_START**

The 1-based starting position of the REFERENCE alignment used in the pairwise alignment between the REFERNECE and QUERY. 

This is derived from alignments of the isoforms earlier against a reference assembly. (Note: This "reference" is unrelated to the REFERENCE used in the pairwise alignment.)

**3. POS**

The 1-based position of the variation on the given sequence.

Same as the VCF format. https://en.wikipedia.org/wiki/Variant_Call_Format

**4. ID**

Pair of IDs used for REFERENCE and QUERY used in the pairwise alignment.

Users will need to investigate these, e.g. to figure which which isoform was used in these pairwise alignments.


**5. REF**

The reference base (or bases in the case of an InDel) at the given position on the given REFERENCE sequence (whereby the "REFERENCE" is the REFERENCE in the pairwise alignment).

Same as the VCF format. https://en.wikipedia.org/wiki/Variant_Call_Format

**6. ALT**

The list of alternative alleles at this position.

Same as the VCF format. https://en.wikipedia.org/wiki/Variant_Call_Format


**7. REFERENCE**

The sequence of the "REFERENCE" when we do the pairwise alignment between the REFERENCE and the QUERY.

**8. QUERY**

The sequence of the "QUERY" when we do the pairwise alignment between the REFERENCE and the QUERY.




