# call_variants.py

## Help message
```
Description: The program calls SNPs and InDels from the unique isoform file

positional arguments:
  inFile       string	input file.  e.g. /in/unique_isoforms.tsv
  outFile      string	output file.  e.g. /out/unique_isoforms_variants.tsv

optional arguments:
  -h, --help   show this help message and exit
  -n , --N     integer	maximum bases for indel combination (inclusive), default 10
  -d, --debug  boolean	show debug information.
```

## Running script
```
python3 call_variants.py call_variants_input.tsv call_variants_output.tsv
```

## Output
CHR: chromosome (contig) that the reference isoform is aligned to
REF_START: starting position of chromosome (contig) that the reference isoform is aligned to
POS: variant starting position on reference isoform sequence
ID: alignment ID in "ref_isoform-query_isoform" format
REF: variant sequence in reference isoform
ALT: altered sequence in query isoform
REFERENCE: full reference isoform sequence
QUERY: full query isoform sequence
```
##This is a header
#CHR	REF_START	POS	ID	REF	ALT	REFERENCE	QUERY
NC_060925.1	255185	4	HG002_PB.6.2-HG004_PB.6.2	G	TTATCCGGAG	GGAGCCAAGGTCCGCTCGGGTGAGTGCCCTCCGCTTTTTGTGGCCAAACCCAGCCACGCAGTCCCCTCCCTGCGGCGTCCTCCACACCCGGGGTCTGCTGGTCTCCGCGGATGTCACAGGCTCGGCAACCGCCCTCCTGTCGGCGGGGAGTCCCGCGACGCCCGGAAATGCTCCGAAGCCTGTCGCTCAGCTGCCAGATCTGCGTCTGTGTCCGGTTCCGTCACTGAGGTCGCCCCTGTCCGGCCCTTCCACCCTAGTTCTCTTCACCGTCCGCCCATCCTATCGCGCGCGGCCTCAGGTCCCCATTCGGCATGTGGCTTGTCTTCCATCGTCCCCACCCTCGCCCCTCTTGGCCCCTCAGGGCAGCCCTGGGATTCGGCAGACGCCAGTCCTCCCTGAGATGCTTCCCCGTCCTTCCCTCCGCCAGGCCCTACGTCTCCGCAAACCCCACGCTTCGGGGTGGCCGCCTCAGACAGGACCCTGAGTCCGAGACTGGGGTAGGGGACCTGCCCGATCCTGTAACAACCCTCGTGCTTCTGCACAATCGCCTCCCACTAGCGGTGACTGTTGGGTGTCTACCTTCCCGGTGTCCCACTGAGAAGCGGGCTCCTCCTTGGCAGGGGCTTCTTCATTGCCTCGCTGTGGATGTCGAGGTGGGGCAGGAGAGTGAGGAGAAAACAGAGAGGAGGGAGGTAGAGCCAACGAGCGAGAAAAGGGGAGGGAAGTTTAGATGGGAAGTGGATGGGTCTGAGGAATTTGAACAAACACCGACAGTGAAGGAGAGTGACCTGAGCAAGCAGTAGTGGGGTAAATGGAAATAGACAAAATGGGAATCAGCAGAGATATGGAGGACAGAATACAATGAGGAGGCCTTGACCGTCAGTAGCAGAGAGGGCAGCAGAAGCCTAATTCCCAAATTCCTTAGTGGTTTTCTGATTTCCAAATTAGTTTCCCTTTTAAATTTATTGTGTCAGGTTCAGCTTATGAGGCCTCAATACTTTTCAGTCTTAATTGTATATTGAAAATACTTTTTGTTTACTAAATGCTTTTTACATTAATTCAGTGTGCACTCCGTAAGGATATTGATGATTTGAGTTAGTTTAGTATTCAACAGCTTCCTCTATTCCTTTGTATGATCTCTGTATTTAATGGCTGTGGCATAAAGTTTCCAACTAAGTATAAGTATCAAGTTTTCTTTGTGATGTTTTCTGCAAATATTGAAGGATGACCTGGATTGTCCTAGAACTTTGTTCCAACAGATTACATGTGTTCATAACGAATAAATTGCTC	GGATTATCCGGAGCCAAGGTCCGCTCGGGTGAGTGCCCTCCGCTTTTTGTGGCCAAACCCAGCCACGCAGTCCCCTCCCTGCGGCGTCCTCCACACCCGGGGTCTGCTGGTCTCCGCGGATGTCACAGGCTCGGCAACCGCCCTCCTGTCGGCGGGGAGTCCCGCGACGCCCGGAAATGCTCCGAAGCCTGTCGCTCAGCTGCCAGATCTGCGTCTGTGTCCGGTTCCGTCACTGAGGTCGCCCCTGTCCGGCCCTTCCACCCTAGTTCTCTTCACCGTCCGCCCATCCTATCGCGCGCGGCCTCAGGTCCCCATTCGGCATGTGGCTTGTCTTCCATCGTCCCCACCCTCGCCCCTCTTGGCCCCTCAGGGCAGCCCTGGGATTCGGCAGACGCCAGTCCTCCCTGAGATGCTTCCCCGTCCTTCCCTCCGCCAGGCCCTACGTCTCCGCAAACCCCACGCTTCGGGGTGGCCGCCTCAGACAGGACCCTGAGTCCGAGACTGGGGTAGGGGACCTGCCCGATCCTGTAACAACCCTCGTGCTTCTGCACAATCGCCTCCCACTAGCGGTGACTGTTGGGTGTCTACCTTCCCGGTGTCCCACTGAGAAGCGGGCTCCTCCTTGGCAGGGGCTTCTTCATTGCCTCGCTGTGGATGTCGAGGTGGGGCAGGAGAGTGAGGAGAAAACAGAGAGGAGGGAGGTAGAGCCAACGAGCGAGAAAAGGGGAGGGAAGTTTAGATGGGAAGTGGATGGGTCTGAGGAATTTGAACAAACACCGACAGTGAAGGAGAGTGACCTGAGCAAGCAGTAGTGGGGTAAATGGAAATAGACAAAATGGGAATCAGCAGAGATATGGAGGACAGAATACAATGAGGAGGCCTTGACCGTCAGTAGCAGAGAGGGCAGCAGAAGCCTAATTCCCAAATTCCTTAGTGGTTTTCTGATTTCCAAATTAGTTTCCCTTTTAAATTTATTGTGTCAGGTTCAGCTTATGAGGCCTCAATACTTTTCAGTCTTAATTGTATATTGAAAATACTTTTTGTTTACTAAATGCTTTTTACATTAATTCAGTGTGCACTCCGTAAGGATATTGATGATTTGAGTTAGTTTAGTATTCAACAGCTTCCTCTATTCCTTTGTATGATCTCTGTATTTAATGGCTGTGGCATAAAGTTTCCAACTAAGTATAAGTATCAAGTTTTCTTTGTGATGTTTTCTGCAAATATTGAAGGATGACCTGGATTGTCCTAGAACTTTGTTCCAACAGATTACATGTGTTCATAACGAATAAATTGCTCAAAGATATTTC
NC_060925.1	255185	1291	HG002_PB.6.2-HG004_PB.6.2	.	AAAGATATTTC	GGAGCCAAGGTCCGCTCGGGTGAGTGCCCTCCGCTTTTTGTGGCCAAACCCAGCCACGCAGTCCCCTCCCTGCGGCGTCCTCCACACCCGGGGTCTGCTGGTCTCCGCGGATGTCACAGGCTCGGCAACCGCCCTCCTGTCGGCGGGGAGTCCCGCGACGCCCGGAAATGCTCCGAAGCCTGTCGCTCAGCTGCCAGATCTGCGTCTGTGTCCGGTTCCGTCACTGAGGTCGCCCCTGTCCGGCCCTTCCACCCTAGTTCTCTTCACCGTCCGCCCATCCTATCGCGCGCGGCCTCAGGTCCCCATTCGGCATGTGGCTTGTCTTCCATCGTCCCCACCCTCGCCCCTCTTGGCCCCTCAGGGCAGCCCTGGGATTCGGCAGACGCCAGTCCTCCCTGAGATGCTTCCCCGTCCTTCCCTCCGCCAGGCCCTACGTCTCCGCAAACCCCACGCTTCGGGGTGGCCGCCTCAGACAGGACCCTGAGTCCGAGACTGGGGTAGGGGACCTGCCCGATCCTGTAACAACCCTCGTGCTTCTGCACAATCGCCTCCCACTAGCGGTGACTGTTGGGTGTCTACCTTCCCGGTGTCCCACTGAGAAGCGGGCTCCTCCTTGGCAGGGGCTTCTTCATTGCCTCGCTGTGGATGTCGAGGTGGGGCAGGAGAGTGAGGAGAAAACAGAGAGGAGGGAGGTAGAGCCAACGAGCGAGAAAAGGGGAGGGAAGTTTAGATGGGAAGTGGATGGGTCTGAGGAATTTGAACAAACACCGACAGTGAAGGAGAGTGACCTGAGCAAGCAGTAGTGGGGTAAATGGAAATAGACAAAATGGGAATCAGCAGAGATATGGAGGACAGAATACAATGAGGAGGCCTTGACCGTCAGTAGCAGAGAGGGCAGCAGAAGCCTAATTCCCAAATTCCTTAGTGGTTTTCTGATTTCCAAATTAGTTTCCCTTTTAAATTTATTGTGTCAGGTTCAGCTTATGAGGCCTCAATACTTTTCAGTCTTAATTGTATATTGAAAATACTTTTTGTTTACTAAATGCTTTTTACATTAATTCAGTGTGCACTCCGTAAGGATATTGATGATTTGAGTTAGTTTAGTATTCAACAGCTTCCTCTATTCCTTTGTATGATCTCTGTATTTAATGGCTGTGGCATAAAGTTTCCAACTAAGTATAAGTATCAAGTTTTCTTTGTGATGTTTTCTGCAAATATTGAAGGATGACCTGGATTGTCCTAGAACTTTGTTCCAACAGATTACATGTGTTCATAACGAATAAATTGCTC	GGATTATCCGGAGCCAAGGTCCGCTCGGGTGAGTGCCCTCCGCTTTTTGTGGCCAAACCCAGCCACGCAGTCCCCTCCCTGCGGCGTCCTCCACACCCGGGGTCTGCTGGTCTCCGCGGATGTCACAGGCTCGGCAACCGCCCTCCTGTCGGCGGGGAGTCCCGCGACGCCCGGAAATGCTCCGAAGCCTGTCGCTCAGCTGCCAGATCTGCGTCTGTGTCCGGTTCCGTCACTGAGGTCGCCCCTGTCCGGCCCTTCCACCCTAGTTCTCTTCACCGTCCGCCCATCCTATCGCGCGCGGCCTCAGGTCCCCATTCGGCATGTGGCTTGTCTTCCATCGTCCCCACCCTCGCCCCTCTTGGCCCCTCAGGGCAGCCCTGGGATTCGGCAGACGCCAGTCCTCCCTGAGATGCTTCCCCGTCCTTCCCTCCGCCAGGCCCTACGTCTCCGCAAACCCCACGCTTCGGGGTGGCCGCCTCAGACAGGACCCTGAGTCCGAGACTGGGGTAGGGGACCTGCCCGATCCTGTAACAACCCTCGTGCTTCTGCACAATCGCCTCCCACTAGCGGTGACTGTTGGGTGTCTACCTTCCCGGTGTCCCACTGAGAAGCGGGCTCCTCCTTGGCAGGGGCTTCTTCATTGCCTCGCTGTGGATGTCGAGGTGGGGCAGGAGAGTGAGGAGAAAACAGAGAGGAGGGAGGTAGAGCCAACGAGCGAGAAAAGGGGAGGGAAGTTTAGATGGGAAGTGGATGGGTCTGAGGAATTTGAACAAACACCGACAGTGAAGGAGAGTGACCTGAGCAAGCAGTAGTGGGGTAAATGGAAATAGACAAAATGGGAATCAGCAGAGATATGGAGGACAGAATACAATGAGGAGGCCTTGACCGTCAGTAGCAGAGAGGGCAGCAGAAGCCTAATTCCCAAATTCCTTAGTGGTTTTCTGATTTCCAAATTAGTTTCCCTTTTAAATTTATTGTGTCAGGTTCAGCTTATGAGGCCTCAATACTTTTCAGTCTTAATTGTATATTGAAAATACTTTTTGTTTACTAAATGCTTTTTACATTAATTCAGTGTGCACTCCGTAAGGATATTGATGATTTGAGTTAGTTTAGTATTCAACAGCTTCCTCTATTCCTTTGTATGATCTCTGTATTTAATGGCTGTGGCATAAAGTTTCCAACTAAGTATAAGTATCAAGTTTTCTTTGTGATGTTTTCTGCAAATATTGAAGGATGACCTGGATTGTCCTAGAACTTTGTTCCAACAGATTACATGTGTTCATAACGAATAAATTGCTCAAAGATATTTC
```
