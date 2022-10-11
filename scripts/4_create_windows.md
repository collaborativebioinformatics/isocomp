## Create windows to prepare for sequence extraction from binning

```
USAGE:  python scripts/4_create_windows_yq.py <bp_cov_names> <gene_window_name>

        bp_cov_names            comma-delimited paths to 100bp windows with coverage for each sample
                                   expected fields: [chr, start, end, per base coverage]
        gene_window_name        path to transcript/RNA boundaries
                                   expected fields: [chr, start, end, geneName, 0, -, geneID, biotype]
```

### Input
1. 100bp windows with coverage for each sample
    - expected fields: ["chr", "start", "end", "per base coverage"]
    - example:
    ```
    NC_060925.1     255501  255600  1.00
    NC_060925.1     255601  255700  1.00
    NC_060925.1     255701  255800  1.00
    NC_060925.1     255801  255900  1.00
    NC_060925.1     255901  256000  1.00
    ```


2. transcript/RNA boundaries
    - expected fields: `["chr", "start", "end", "geneName", "0", "-", "geneID", "biotype"]`
    - example:
    ```
    NC_060943.1     61441599        61449907        A1BG    0       -       MIM:138670      protein_coding
    NC_060943.1     61448385        61451599        A1BG-AS1        0       +       HGNC:HGNC:37133 lncRNA
    NC_060934.1     51648044        51734261        A1CF    0       -       MIM:618199      protein_coding
    ```

### Output
#### Intermediate files
(Can be deleted after running the script)
- `intermediate_files/HG_merged_100bp_coverage.bed` 
    - merged 100bp coverage bed files from input 1.
- `intermediate_files/merged_100bp_coverage.nonzero.bed`
    - intersected bed between gene windows (input 2) and 100 bp windows (input 1) along the entire genome

#### Output files
- `output_files/4_merged_windows.bed`
    - **THE ONLY OUTPUT THAT IS USED AS INPUT FOR THE NEXT STEP**
    - merged windows between gene boundaries and 100bp coverage windows
    - columns: ["chr", "start", "end"]
    - example:
    ```
    NC_060925.1	450247	486181
    NC_060925.1	510436	544704
    NC_060925.1	645388	660600
    ```
- `output_files/gene_overlap_sum_coverage.bed`
    - sum of coverage in 100bp bins for each gene.
    - columns: ["chr", "start", "end", "coverage", "geneName"]
    - example:
    ```
    NC_060925.1	255178	256494	25.57	LINC00115
    NC_060925.1	256563	288416	25.2	LINC01128
    NC_060925.1	372945	388041	209.14	NOC2L
    NC_060925.1	353566	373316	31.21	SAMD11
    ```
- `output_files/no_overlap_100bp.nonzero.bed`
    - 100 bp windows with nonzero coverage with no overlap with genes
    - columns: ["chr", "start","end", "hg002", "hg004", "hg005", "sum coverage"]
    - example"
    ```
    NC_060925.1	660501	660600	0.0	0.14	0.1	0.24
    NC_060925.1	1093501	1093600	1.0	1.0	0.57	2.57
    NC_060925.1	1093601	1093700	0.47	0.94	0.0	1.41
    ```