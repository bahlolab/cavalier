#' @importFrom readr read_delim cols
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
      setNames(c('symbol', 'RVIS_percentile')) %>%
      mutate(symbol = if_else(symbol %in% hgnc_alias$symbol,
                              symbol,
                              hgnc_alias$symbol[match(symbol, hgnc_alias$alias)])) %>% 
      filter(!is.na(symbol))
    
    options('cavalier.rvis_table' = rvis_table)
  }
  
  return(rvis_table)
}