
#' Retrieve GeVIR table from get_cavalier_opt("gevir_url") or disk cache
#' @importFrom readr read_csv
#' @importFrom dplyr "%>%" mutate rename select
get_gevir_table <- function()
{
  gevir_url <- get_cavalier_opt("gevir_url")
  
  fun <- function() {
    retry('GET', gevir_url) %>% 
      content(as = 'raw') %>% 
      rawConnection() %>% 
      read_csv(col_types = cols()) %>% 
      select(symbol = gene_name, ensembl_gene_id = gene_id, gevir_percentile, loeuf_percentile) %>% 
      distinct()
  }
  
  cache(
    fun = fun,
    name = 'gevir_gene_rankings'
  )
}

#' Update GeVIR table gene symbol with current HGNC symbol
get_gevir_table_hgnc <- function()
{
  get_gevir_table() %>% 
    mutate(symbol = coalesce(hgnc_ensembl2sym(ensembl_gene_id),
                             hgnc_sym2sym(symbol)))
}

#' Convert HGNC gene symbol to GeVIR percentile
#' @export
sym2gevir <- function(symbol)
{
  symbol <- hgnc_sym2sym(symbol)
  gevir_table <- get_gevir_table_hgnc()
  gevir_table$gevir_percentile[match(symbol, gevir_table$symbol)]
}

#' Convert ensemble gene id to GeVIR percentile
#' @export
ensembl2gevir <- function(ensembl)
{
  gevir_table <- get_gevir_table_hgnc()
  gevir_table$gevir_percentile[match(ensembl, gevir_table$ensembl_gene_id)]
}

#' Convert HGNC gene symbol to (ExAC) LOUEF percentile 
#' @export
sym2loeuf <- function(symbol)
{
  symbol <- hgnc_sym2sym(symbol)
  gevir_table <- get_gevir_table_hgnc()
  gevir_table$loeuf_percentile[match(symbol, gevir_table$symbol)]
}

#' Convert ensemble gene id to to (ExAC) LOUEF percentile 
#' @export
ensembl2loeuf <- function(ensembl)
{
  gevir_table <- get_gevir_table_hgnc()
  gevir_table$loeuf_percentile[match(ensembl, gevir_table$ensembl_gene_id)]
}
