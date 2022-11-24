# Trim indel ends

## Description
This program looks at pair-wise isoform alignments within each window and trims off beginning and ending of isoforms to reduce false posistives.

## Help message
`python scripts/trim_indel_ends.py manifest unique_isoform_aln [-m mode -t threshold]`

```
- manifest -- manifest file in form <sample,in/fasta.fa> that describes the sample and location of fasta files
- unique_isoform_aln -- output from `find_unique_isoform.py`
- mode -- optional -- "strict" or "lenient"
    - "strict" will only trim if after trimming, two isoforms are exactly the same
    - "lenient" will trim whenever indels are encountered at the ends of the alignment
    - Default is "lenient"
- threshold -- optional -- maximum number of bases to be trimmed off each end. Default is 10000 (no limit)
```

## Input
- fasta files for each sample
- `unique_isoform_and_aln_stats.tsv` produced by `find_unique_isoform.py`

## Output
- separate sample fasta files
- Histogram describing trim lengths
- Scatter plot describing edit distance before and after trimmming

## Logic
1. For each isoform present in `unique_isoform_and_aln_stats.tsv`, this program first finds the isoform in another sample that is the most similar to it. 

2. Then, this program looks at the cigar string of each pair of sequences. 
    - If the cigar string starts or ends with insertions `"I"` or deletions `"D"`, the program trims the sequences.
        - Trim the first-in-pair if `"I"` is encountered
        - Trim the second-in-pair if `"D"` is encountered
3. Write sequences in separate sample fasta files
    - If a sequence is trimmed in multiple pairs, keep the longest one
    - If a sequence is not included in `unique_isoform_and_aln_stats.tsv`, look for it in original fasta and include it.

## Dependencies
```
import pandas as pd
from collections import defaultdict, Counter
import re
import edlib
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import matplotlib.pyplot as plt
from tqdm import tqdm
import os
import sys
import argparse
import datetime
```