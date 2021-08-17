
#' @importFrom readr cols read_tsv
#' @importFrom dplyr "%>%" mutate rename
#' @importFrom stringr str_remove
#' @importFrom rlang is_scalar_character
#' @export
get_hgnc_complete <- function()
{
    hgnc_complete <- getOption('cavalier.hgnc_complete')
    
    if (is.null(hgnc_complete)) {
        
        cache_dir <- getOption('cavalier.cache_dir')
        hgnc_complete_uri <- getOption('cavalier.hgnc_complete_uri')
        fn <- basename(hgnc_complete_uri) %>% 
            str_remove('.txt$') %>% 
            str_c('.rds') %>% 
            file.path(cache_dir, .)
        
        hgnc_complete <- 
            suppressWarnings(
                (function() read_tsv(hgnc_complete_uri, col_types = cols())) %>% 
                    cache(fn)) %>% 
            select(symbol, name, location, ensembl_gene_id, alias_symbol, prev_symbol)
        
        options('cavalier.hgnc_complete' = hgnc_complete)
    }
    return(hgnc_complete)
}

#' @importFrom tidyr replace_na separate_rows
#' @importFrom stringr str_c
#' @importFrom dplyr "%>%" mutate rename select if_else add_count filter case_when
#' @export
get_hgnc_alias <- function() 
{
    hgnc_alias <- getOption('cavalier.hgnc_alias')
    
    if (is.null(hgnc_alias)) {
        
        hgnc_complete <- 
            get_hgnc_complete()
        
        hgnc_alias <-
            hgnc_complete %>% 
            # filter(!(is.na(alias_symbol) & is.na(prev_symbol))) %>% 
            mutate(alias = case_when(
                !is.na(alias_symbol) & !is.na(prev_symbol) ~ str_c(alias_symbol, '|', prev_symbol),
                !is.na(alias_symbol)                       ~ alias_symbol,
                !is.na(prev_symbol)                        ~ prev_symbol)) %>% 
            select(symbol, ensembl_gene_id,  alias) %>% 
            separate_rows(alias, sep = '\\|') %>% 
            mutate(alias = if_else(alias %in% symbol, NA_character_, alias)) %>% 
            add_count(alias) %>%
            mutate(alias = if_else(n != 1, NA_character_, alias)) %>% 
            select(-n) %>% 
            distinct() %>% 
            add_count(symbol) %>% 
            filter(!(is.na(alias) & n > 1)) %>% 
            select(-n)
        
        options('cavalier.hgnc_alias' = hgnc_alias)
    }
    
    return(hgnc_alias)
}

# replace gene symbols alias with HGNC approved symbol
#' @export
hgnc_name_replace <- function(genes) {
    hgnc_alias <- get_hgnc_alias()
    at <- which(genes %in% hgnc_alias$alias)
    replace(genes, at, hgnc_alias$symbol[match(genes[at], hgnc_alias$alias)])
}
