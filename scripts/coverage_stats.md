
# coverage_stats.py

## Help message
```
Description: The program generate bam coverage statistics for given bed regions

positional arguments:
  inBam            string	input bam file.  e.g. /in/sample.bam
  inBed            string	input bed file.  e.g. /in/sample.bed
  samtools         string	samtools.  e.g. /path/to/samtools
  outBed           string	output regional coverage statistics.  e.g. /out/sample.bed

optional arguments:
  -h, --help       show this help message and exit
  -q , --mapQCut   integer	minimum mapQ for counted reads, default 20

```

## Running script
```
python3 bamStats.py sorted_hg002.bam T2T.100kb.bed /path/to/samtools HG002.bed -q 0
```

## Output
Columns: window_chr, window_start, window_end, average_depth
```
NC_060925.1     1       100     0.00
NC_060925.1     101     200     0.00
NC_060925.1     201     300     0.00
```

