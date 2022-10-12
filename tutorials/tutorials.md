# Create Genomic Bins

1. Download the reference annotations for your genome. [See here the T2T CHM13v2.0 
genome](https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_009914755.1/)

2. Use [GenomicFeatures package](https://www.bioconductor.org/packages/release/bioc/html/GenomicFeatures.html) to parse the gtf/gff into a sqlite database

```R
txdb = GenomicFeatures::makeTxDbFromGFF('path/to/genome.gtf', format='gtf')
```

3. Another method for parsing the gtf/gff is with [rtracklayer](https://www.bioconductor.org/packages/release/bioc/html/rtracklayer.html). This 
will contain all the information in the annotation file, including that which isn't 
parsed into the txdb automatically. It can be used to manually augment the 
txdb. I use both here, but the gff_granges object alone would be sufficient.

```R
gff_granges = rtracklayer::import('path/to/genome.gtf')
```

4. To extract the gene regions for binning, we use the following:

```R
  library(tidyverse)
  library(GenomicFeatures) # bioconductor
  library(regioneR) # bioconductor
  
  date = format(Sys.Date(), format="%Y%m%d")
  
  # Extract all gene features to a bed6+
  all_gene_regions = genes(txdb) %>%
    as_tibble() %>%
    left_join(
      as_tibble(gff_granges) %>%
        dplyr::select(gene_id, db_xref, gene, gene_biotype) %>%
        distinct(gene_id, .keep_all = TRUE),
      by = 'gene_id') %>%
    # to comply with bed6 format
    mutate(score = 0) %>%
    select(seqnames,start,end,gene,score,strand,db_xref,gene_biotype) %>%
    write_tsv(sprintf("all_gene_regions_%s.bed",date),
              col_names = FALSE)

  # Extract only the protein coding regions to a bed6+ file
  protein_coding_gene_regions = genes(txdb) %>%
    as_tibble() %>%
    left_join(
      as_tibble(gff_granges) %>%
        dplyr::select(gene_id, db_xref, gene, gene_biotype) %>%
        distinct(gene_id, .keep_all = TRUE),
      by = 'gene_id') %>%
    # to comply with bed6 format
    mutate(score = 0) %>%
    select(seqnames,start,end,gene,score,strand,db_xref,gene_biotype) %>%
    filter(gene_biotype == 'protein_coding') %>%
    write_tsv(sprintf("protein_coding_regions_%s.bed",date),
              col_names = FALSE)

  # intersect overlapping protein coding features, export in bed3 format
  collapsed_protein_coding_gene_regions = protein_coding_gene_regions %>%
    dplyr::rename(chr=seqnames) %>%
    as.data.frame() %>%
    regioneR::toGRanges() %>%
    regioneR::joinRegions(min.dist = 1) %>%
    regioneR::toDataframe() %>%
    as_tibble() %>%
    write_tsv(sprintf("collapsed_protein_coding_%s.bed",date),
              col_names = FALSE)
```
