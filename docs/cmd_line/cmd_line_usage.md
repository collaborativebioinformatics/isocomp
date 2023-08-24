# Interface

After installing `isocomp`, you may use the command line interface like so:

```bash
isocomp --help
```
which will produce a help mesage like so:

```raw
usage: isocomp [-h] [-v] {create_windows,find_unique_isoforms} ...

isocomp: 0.1.0

positional arguments:
  {create_windows,find_unique_isoforms}
                        Available Tools
    create_windows      from a list of SQANTI gtfs, create a bed6 formatted file -- 0 indexed, half open intervals -- where each entry in the bed file represents a discrete region of
                        coverage over the reference genome.
    find_unique_isoforms
                        Using the clustered gtf output of create_windows and a csv file with columns 'source' and 'fasta' which provide a map from the unique values in 'Source' in the
                        clustered gtf to the corresponding fasta file, find unique isoforms in the clusters described in the clustered_gtf.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
```

Similarly, you can get help for any available tool, for instance

```bash
isocomp create_windows --help
```

produces

```raw
usage: create_windows [-h] [-l {critical,error,warning,info,debug}] -i INPUT [INPUT ...] [-p OUTPUT_PREFIX] [--overwrite | --no-overwrite]

from a list of SQANTI gtfs, create a bed6 formatted file -- 0 indexed, half open intervals -- where each entry in the bed file represents a discrete region of coverage over the
reference genome.

optional arguments:
  -h, --help            show this help message and exit

general:
  -l {critical,error,warning,info,debug}, --log_level {critical,error,warning,info,debug}

input:
  -i INPUT [INPUT ...], --input INPUT [INPUT ...]
                        a space delimited list of paths to gtf files from a concatenated merged set of regions will be created

output:
  -p OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        This should be either a file basename -- NO extension -- which will be written to the current working directory, or a valid path up to a filename, eg
                        /output/path/concat would result in the file output/path/concat.gtf being written
  --overwrite, --no-overwrite
                        By default, create_windows will overwrite a file with the same name in the output directory. Set --no-overwrite to avoid overwriting.
```

# create_windows

If you have a directory which looks like the following:

```raw
├── clustered_regions.gtf
├── fasta_map.csv
├── hg002_classification.txt
├── hg002_sqanti_fltr.fasta
├── hg002_sqanti_fltr.fasta.fai
├── hg002_sqanti_fltr.gtf
├── hg004_classification.txt
├── hg004_sqanti_fltr.fasta
├── hg004_sqanti_fltr.fasta.fai
├── hg004_sqanti_fltr.gtf
├── hg005_classification.txt
├── hg005_sqanti_fltr.fasta
├── hg005_sqanti_fltr.fasta.fai
├── hg005_sqanti_fltr.gtf
└── unique_isoforms.csv
```

Note that this happens to be the test data directory in `isocomp`

Then you could run `create_windows` like so:

```bash
isocomp create_indows -i /path/to/input_dir/*gtf
```

or, more explicitely

```bash
isocomp create_indows -i /path/to/input_dir/hg002_sqanti_fltr.gtf /path/to/input_dir/hg004_sqanti_fltr.gtf /path/to/input_dir/hg005_sqanti_fltr.gtf
```

This will output clustered_regions.gtf in the `$PWD`

# find_unique_isoforms

Input to find_unique_isoforms is clustered gtf which is output by 
`create_windows` and a `csv` file with the columns `source` and `fasta` 
where the records store a value which corresponds to a unique factor level 
in the `clustered_regions.gtf` `Source` column and a path to a fasta file 
which stores the isoform sequences.  

For example, the clustered_gtf might look like this:

```raw
➜  data git:(oop_ify) ✗ head clustered_regions.gtf 
chr1	hg004_sqanti_fltr	transcript	1013497	1014531	.	+	.transcript_id "PB.13.1"; gene_id "PB.13"; Cluster "1";
chr1	hg005_sqanti_fltr	transcript	1013497	1014531	.	+	.transcript_id "PB.17.1"; gene_id "PB.17"; Cluster "1";
chr1	hg005_sqanti_fltr	transcript	1013497	1014531	.	+	.transcript_id "PB.17.2"; gene_id "PB.17"; Cluster "1";
chr1	hg002_sqanti_fltr	transcript	1013504	1014531	.	+	.transcript_id "PB.17.2"; gene_id "PB.17"; Cluster "1";
chr1	hg004_sqanti_fltr	transcript	1013532	1014531	.	+	.transcript_id "PB.13.2"; gene_id "PB.13"; Cluster "1";
chr1	hg005_sqanti_fltr	transcript	1020120	1056112	.	+	.transcript_id "PB.18.3"; gene_id "PB.18"; Cluster "2";

```

while the fasta map would look like this:

```raw
source,fasta
hg002_sqanti_fltr,/home/oguzkhan/code/isocomp/src/tests/data/hg002_sqanti_fltr.fasta
hg004_sqanti_fltr,/home/oguzkhan/code/isocomp/src/tests/data/hg004_sqanti_fltr.fasta
hg005_sqanti_fltr,/home/oguzkhan/code/isocomp/src/tests/data/hg005_sqanti_fltr.fasta

```

To find the unique isoforms in each cluster, use `find_unique_isoforms` as follows:

```bash
isocomp find_unique_isoforms -a clustered_regions.gtf -f fasta_map.csv
```