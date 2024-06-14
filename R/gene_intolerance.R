
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
      mutate(
        symbol = coalesce(hgnc_ensembl2sym(ensembl_gene_id),
                          hgnc_sym2sym(symbol))
      ) %>% 
      distinct()
  }
  
  cache(
    fun = fun,
    name = basename(gevir_url),
    disk = TRUE
  )
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
