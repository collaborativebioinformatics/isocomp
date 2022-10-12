(Shall we use this readme as a list of requirements?)

# create_tx_start_stop_bed.R
The R script for creating tx regions from the T2T chm13v2.0 gtf

## Dependencies
- optparse: cmd line argument parsing
- GenomicFeatures: handling range/feature data eg from a gtf
- dplyr: dataframe handling
- readr: read/write dataframes
- regioneR: functions on region sets, eg union over overlapping regions

# Create windows for binning isoforms

## Dependencies
- python3.9
- pybedtools
- bedtools
- pandas