
# call_variants.py

## Help message
```
Description: The program calls SNPs and InDels from the unique isoform file

positional arguments:
  inFile      string	input file.  e.g. /in/unique_isoforms.tsv
  outFile     string	output file.  e.g. /out/unique_isoforms_variants.tsv

optional arguments:
  -h, --help  show this help message and exit
  -n , --N    float	maximum bases for indel combination (inclusive), default 10
```

## Running script
```
python3 call_variants.py unique_isoforms_and_aln_stats_filter0.99.tsv unique_isoforms_and_aln_stats_filter0.99.call.tsv --N 10
```

## Output
ref_id: reference sequence ID
query_id: query sequence ID
type: type of variant (snp/ins/del/delins)
ref: variant sequence in reference
alt: altered sequence
ref_start: starting position of the variant in reference
```
#ref_id	query_id	type	ref	alt	ref_pos
HG002_PB.6.2	HG004_PB.6.2	delins	G	GCCAAGGTCC	4
HG002_PB.6.2	HG004_PB.6.2	ins	.	AAAGATATTTC	1291
HG002_PB.6.2	HG005_PB.7566.3	snp	A	G	7
```

