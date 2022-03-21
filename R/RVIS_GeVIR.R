
rvis_uri <- "http://genic-intolerance.org/data/RVIS_Unpublished_ExACv2_March2017.txt"

#' @importFrom readr read_tsv cols
#' @importFrom dplyr "%>%" mutate rename
get_rvis_table <- function()
{
  (function()
    retry('GET', rvis_uri) %>% 
     content(as = 'raw') %>% 
     rawConnection() %>% 
     read_tsv(col_types = cols()) %>% 
     select(1,4) %>% 
     setNames(c('symbol', 'rvis_percentile')) %>%
     mutate(symbol = hgnc_sym2sym(symbol))) %>% 
    cache(str_remove(basename(rvis_uri), '.txt$'),
          disk = TRUE)
}

sym2rvis <- function(symbol)
{
  rvis_table <- get_rvis_table()
  rvis_table$rvis_percentile[match(symbol, rvis_table$symbol)]
}

ensembl2rvis <- function(ensembl)
{
  rvis_table <- get_rvis_table()
  rvis_table$rvis_percentile[match(hgnc_ensembl2sym(ensembl), rvis_table$symbol)]
}

entrez2rvis <- function(entrez)
{
  rvis_table <- get_rvis_table()
  rvis_table$rvis_percentile[match(hgnc_entrez2sym(entrez), rvis_table$symbol)]
}

gevir_uri <- "http://www.gevirank.org/static/files/gene_ranking.csv"
#' @importFrom readr read_csv
#' @importFrom dplyr "%>%" mutate rename select
get_gevir_table <- function()
{
  (function() 
    retry('GET', gevir_uri) %>% 
     content(as = 'raw') %>% 
     rawConnection() %>% 
     read_csv(col_types = cols()) %>% 
     mutate(symbol = coalesce(hgnc_ensembl2sym(gene_id),
                              hgnc_sym2sym(gene_name))) %>% 
     select(symbol, ensembl_gene_id = gene_id, gevir_percentile, loeuf_percentile)) %>% 
    cache(str_remove(basename(gevir_uri), '.csv$'),
          disk = TRUE)
}

sym2gevir <- function(symbol)
{
  gevir_table <- get_gevir_table()
  gevir_table$gevir_percentile[match(symbol, gevir_table$symbol)]
}

ensembl2gevir <- function(ensembl)
{
  gevir_table <- get_gevir_table()
  gevir_table$gevir_percentile[match(ensembl, gevir_table$ensembl_gene_id)]
}

entrez2gevir <- function(entrez)
{
  gevir_table <- get_gevir_table()
  gevir_table$gevir_percentile[match(hgnc_entrez2sym(entrez), gevir_table$symbol)]
}

sym2loeuf <- function(symbol)
{
  gevir_table <- get_gevir_table()
  gevir_table$loeuf_percentile[match(symbol, gevir_table$symbol)]
}

ensembl2loeuf <- function(ensembl)
{
  gevir_table <- get_gevir_table()
  gevir_table$loeuf_percentile[match(ensembl, gevir_table$ensembl_gene_id)]
}

entrez2loeuf <- function(entrez)
{
  gevir_table <- get_gevir_table()
  gevir_table$loeuf_percentile[match(hgnc_entrez2sym(entrez), gevir_table$symbol)]
}

gnomad_constraint_uri <- "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"
#' @importFrom readr read_csv
#' @importFrom dplyr "%>%" mutate rename select
get_gnomad_constraints <- function()
{
  (function() {
    tmp <- tempfile(fileext = '.gz')
    download.file(gnomad_constraint_uri, tmp)
    x <- 
      read_tsv(tmp, col_types = cols()) %>% 
      mutate(symbol = hgnc_ensembl2sym(gene_id)) %>% 
      select(symbol, gene_id, transcript, oe_lof, oe_lof_lower, oe_lof_upper, oe_mis, oe_mis_lower, oe_mis_upper, pLI)
    file.remove(tmp) 
    return(x)
    }) %>% 
    cache(str_remove(basename(gnomad_constraint_uri), '.txt.bgz$'),
          disk = TRUE)
}

sym2oe_lof <- function(symbol)
{
  gnomad_constraints <- get_gnomad_constraints()
  gnomad_constraints$oe_lof[match(symbol, gnomad_constraints$symbol)]
}

ensembl2oe_lof <- function(ensembl)
{
  gnomad_constraints <- get_gnomad_constraints()
  gnomad_constraints$oe_lof[match(ensembl, gnomad_constraints$gene_id)]
}

