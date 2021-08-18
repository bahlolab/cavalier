#' @importFrom readr read_tsv cols
#' @importFrom dplyr "%>%" mutate rename
get_rvis_table <- function()
{
  rvis_table <- getOption('cavalier.rvis_table')
  
  if (is.null(rvis_table)) {
    
    cache_dir <- getOption('cavalier.cache_dir')
    rvis_uri <- getOption('cavalier.rvis_uri')
    fn <- basename(rvis_uri) %>% 
      str_remove('.txt$') %>% 
      str_c('.rds') %>% 
      file.path(cache_dir, .)
    
    hgnc_alias <- get_hgnc_alias()
    
    rvis_table <- 
      (function() read_tsv(rvis_uri, col_types = cols())) %>% 
      cache(fn) %>% 
      select(1,4) %>% 
      setNames(c('symbol', 'rvis_percentile')) %>%
      mutate(symbol = hgnc_sym2sym(symbol))
    
    options('cavalier.rvis_table' = rvis_table)
  }
  
  return(rvis_table)
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

#' @importFrom readr read_csv
#' @importFrom dplyr "%>%" mutate rename select
get_gevir_table <- function()
{
  gevir_table <- getOption('cavalier.gevir_table')
  
  if (is.null(gevir_table)) {
    
    cache_dir <- getOption('cavalier.cache_dir')
    gevir_uri <- getOption('cavalier.gevir_uri')
    fn <- basename(gevir_uri) %>% 
      str_remove('.txt$') %>% 
      str_c('.rds') %>% 
      file.path(cache_dir, .)
    
    gevir_table <- 
      (function() read_csv(gevir_uri, col_types = cols())) %>% 
      cache(fn) %>% 
      mutate(symbol = hgnc_ensembl2sym(gene_id),
             symbol = if_else(is.na(symbol), gene_name, symbol)) %>% 
      select(symbol, ensembl_gene_id = gene_id, gevir_percentile, loeuf_percentile)
    
    options('cavalier.gevir_table' = gevir_table)
  }
  
  return(gevir_table)
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
