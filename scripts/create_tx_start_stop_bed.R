#!/usr/bin/env Rscript

#' Parse a gtf file into objects. Output a txdb sqlite database, a bed6+ file
#' of all gene features, a bed6+ with only protein coding genes, and a bed3
#' with overlapping features combined
#' https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_009914755.1/
#' usage: ./create_tx_start_stop_bed.R -g GCF_009914755.1/genomic.gtf -t txdb.sqlite

suppressMessages(library('optparse'))
suppressMessages(library('RSQLite'))
suppressMessages(library('GenomicFeatures'))
suppressMessages(library('tidyverse'))
suppressMessages(library('regioneR'))

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
    message('constructing txdb...')
    t2t_txdb = makeTxDbFromGFF(path_to_gtf, format='gtf')
    saveDb(t2t_txdb, path_to_txdb)
  # else, load from file
  } else{
    message('connecting to txdb...')
    t2t_txdb = loadDb(path_to_txdb)
  }

  # read in the gtf file as a GRanges obj
  message('loading gtf as GRanges...')
  t2t_gff = rtracklayer::import(path_to_gtf)

  # return the txdb connection and the gff granges obj
  list(
    txdb = t2t_txdb,
    gff_granges = t2t_gff
  )

}

parseArguments <- function() {
  # parse and return cmd line input

  option_list <- list(
    make_option(c('-g', '--gtf'),
                help='path to the gtf annotation set. Suggested source:  https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_009914755.1/'),
    make_option(c('-t', '--txdb'),
                 help='path to a bioconductor sqlite TxDb. If one DNE, it will be created. In that case, provide a writeable path.'))

  args <- parse_args(OptionParser(option_list=option_list))
  return(args)
}

main = function(args){

  date = format(Sys.Date(), format="%Y%m%d")

  if(!file.exists(args$gtf)){
    stop(sprintf('GTF file: %s DNE!', args$gtf))
  }

  if(!file.exists(args$txdb)){
    message(sprintf('txdb: %s DNE. One will be constructed', args$txdb))
    construct_txdb = TRUE
  } else{
    construct_txdb = FALSE
  }

  # parse the annotation resources
  annote_resources =
    construct_annote_resources(args$txdb, args$gtf, construct_txdb)

  message('extracting all genes as bed6+')
  all_gene_regions = genes(annote_resources$txdb) %>%
    as_tibble() %>%
    left_join(
      as_tibble(annote_resources$gff_granges) %>%
        dplyr::select(gene_id, db_xref, gene, gene_biotype) %>%
        distinct(gene_id, .keep_all = TRUE),
      by = 'gene_id') %>%
    # to comply with bed6 format
    mutate(score = 0) %>%
    select(seqnames,start,end,gene,score,strand,db_xref,gene_biotype) %>%
    write_tsv(sprintf("all_gene_regions_%s.bed",date),
              col_names = FALSE)

  message('creating protein coding only bed6+...')
  protein_coding_gene_regions = genes(annote_resources$txdb) %>%
    as_tibble() %>%
    left_join(
      as_tibble(annote_resources$gff_granges) %>%
        dplyr::select(gene_id, db_xref, gene, gene_biotype) %>%
        distinct(gene_id, .keep_all = TRUE),
      by = 'gene_id') %>%
    # to comply with bed6 format
    mutate(score = 0) %>%
    select(seqnames,start,end,gene,score,strand,db_xref,gene_biotype) %>%
    filter(gene_biotype == 'protein_coding') %>%
    write_tsv(sprintf("protein_coding_regions_%s.bed",date),
              col_names = FALSE)

  message('creating union of overlap bed3 from protein coding only set...')
  collapsed_protein_coding_gene_regions = protein_coding_gene_regions %>%
    dplyr::rename(chr=seqnames) %>%
    as.data.frame() %>%
    regioneR::toGRanges() %>%
    regioneR::joinRegions(min.dist = 1) %>%
    regioneR::toDataframe() %>%
    as_tibble() %>%
    write_tsv(sprintf("collapsed_protein_coding_%s.bed",date),
              col_names = FALSE)

  message("Done!")
}

main(parseArguments()) # call main method

# for testing
# input_list = list()
# input_list['gtf'] = 'data/GCF_009914755.1/genomic.sorted.gtf'
# input_list['txdb'] = 'data/t2t_txdb.sqlite'
#
# main(input_list)




