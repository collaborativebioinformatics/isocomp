
# Motivation

The purpose of this document is to describe the diversity of isoforms with 
otherwise exact start/end points

## File specification, same isoforms found

### Columns

**1. ref_sample**

Name of the reference sample. This must be known based on how the bedtools 
intersect file was created

**2. ref_isoform**

Name of the reference isoform. This is the fourth column from the bedtools 
intersect output

**3. query_sample**

Name of the query sample. This is the seventh column of the bedtols intersect 
output

**4. query_isoform**

Name of the query isoform. This is column 11

**4. tx_genomic_length**

end-start

**5. edit_distance**

edit distance between query and ref isoform strings

**6. normalized_edit_distance**

`edit_distance / tx_genomic_length`

**7. cigar**

The cigar string of the alignment between the query and the ref