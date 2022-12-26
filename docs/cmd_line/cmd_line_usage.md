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