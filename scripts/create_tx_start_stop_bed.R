#'
#' Parse a gtf file into objects for easy feature extraction
#'

library(RSQLite)
library(GenomicFeatures)
library(tidyverse)
library(regioneR)

#'
#' create a txdb and granges set from a gtf file
#' @param path_to_txdb path to storage location of txdb -- if not already created,
#' then this is where the txdb will be saved
#' @param path_to_gtf path to the gtf file from which to construct the txdb and
#' granges obj
#' @param construct_txdb Boolean, default to FALSE. Set to TRUE if a txdb of
#' this gtf does not already exist
#'
#' @importfrom AnnotationDb saveDb loadDb
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @importFrom rtracklayer import
#'
#' @return a list of two items, the txdb connection and the gff objects (in memory)
#'
construct_annote_resources = function(path_to_txdb, path_to_gtf, construct_txdb = FALSE){

  # create the txdb if construct_txdb is TRUE
  if(construct_txdb){
    t2t_txdb = makeTxDbFromGFF(path_to_gtf, format='gtf')
    saveDb(t2t_txdb, path_to_txdb)
  # else, load from file
  } else{
    t2t_txdb = loadDb(path_to_txdb)
  }

  # read in the gtf file as a GRanges obj
  t2t_gff = rtracklayer::import(path_to_gtf)

  # return the txdb connection and the gff granges obj
  list(
    txdb = t2t_txdb,
    gff_granges = t2t_gff
  )

}

#'
#' Create a table decribing the variety of transcript biotypes listed in the gtf
#'
#' @param gff_granges granges obj constructed from the gtf
#'
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
#'
#' @return a tibble with fields tx_biotype and count
#'
tx_by_biotype = function(gff_granges){
  gff_granges[gff_granges$type == 'transcript']$transcript_biotype %>%
    factor() %>%
    summary() %>%
    as_tibble(rownames='tx_biotype') %>%
    dplyr::rename(count = value)
}

#'
#' extract tx ranges into a bed6 format table
#'
#' @param txdb a txdb connection
#' @param gff_granges gtf/gff represented as a granges obj
#' @param tx_biotype a biotype corresponding to the table produced by tx_by_biotype
#'
#' @importFrom AnnotationDbi select
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate select filter
#'
#' @return
#'
extract_tx_regions = function(txdb, gff_granges, tx_biotype){

  # extract a set of tx names from the gff
  selected_tx = gff_granges[gff_granges$type == 'transcript' &
                              gff_granges$transcript_biotype == tx_biotype]$transcript_id %>%
    unique()

  n_tx = as.character(length(selected_tx))
  message(sprintf("There are %s unique tx selected by tx_biotype: %s",n_tx,tx_biotype))

  # extract the ranges as a bed file from the txdb
  AnnotationDbi::select(txdb,
                        keys = selected_tx,
                        columns = c('TXCHROM', 'TXSTART', 'TXEND',
                                    'TXNAME', 'TXSTRAND', 'GENEID'),
                        keytype = 'TXNAME') %>%
    as_tibble() %>%
    # add score column to comply with bed6 format
    mutate(score = 0) %>%
    # select in order of bed6 format
    dplyr::select(TXCHROM,TXSTART,TXEND,TXNAME,score,TXSTRAND,GENEID) %>%
    # if there are any NA values, filter them. In the t2t annotations, there
    # are three
    filter(!is.na(TXSTART))
}

#'
#' given a table from extract_tx_regions, collapse overlapping regions
#'
#' @param tx_defns the table output by extract_tx_regions
#' @param min_dist see regioneR::joinRegions. Any pair of regions closer than
#' min.dist bases will be fused in a larger region.
#'
#' @importFrom dplyr rename
#' @importFrom tibble as_tibble
#' @importFrom regioneR toGranges joinRegions toDataFrame
#'
#' @return a dataframe with fields chr start stop of fused regions
#'
collapse_overlapping_regions = function(tx_defns,min_dist=1){
  tx_regions = tx_defns %>%
    dplyr::rename(chr=TXCHROM,start=TXSTART,end=TXEND) %>%
    as.data.frame() %>%
    regioneR::toGRanges() %>%
    regioneR::joinRegions(min.dist = min_dist) %>%
    regioneR::toDataframe() %>%
    as_tibble()
}

# parse the annotation resources
annote_resources =
  construct_annote_resources("data/t2t_txdb.sqlite",
                             "data/GCF_009914755.1/genomic.sorted.gtf")

# create a table describing the diversity and number of tx_biotype in gff
tx_biotype_df = tx_by_biotype(annote_resources$gff_granges)

# extract mRNA tx regions
tx_regions =
  extract_tx_regions(annote_resources$txdb, annote_resources$gff_granges, 'mRNA')

# collapse regions
collapsed_tx_regions = collapse_overlapping_regions(tx_regions)

# write out
write_tsv(tx_regions,
          "data/mrna_tx_start_stop_t2t_20221010.bed",
          col_names = FALSE)

write_tsv(collapsed_tx_regions,
          "data/tx_mrna_nonoverlapping_20221010.bed",
          col_names = FALSE)





