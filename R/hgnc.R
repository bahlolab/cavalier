
#' @importFrom readr cols read_tsv
#' @importFrom dplyr "%>%" mutate rename
#' @importFrom stringr str_remove
#' @importFrom rlang is_scalar_character
get_hgnc_complete <- function()
{
    hgnc_complete <- cavalier_cache$hgnc_complete
    
    if (is.null(hgnc_complete)) {
        
        cache_dir <- get_cavalier_opt('cache_dir')
        hgnc_complete_uri <- get_cavalier_opt('hgnc_complete_uri')
        fn <- basename(hgnc_complete_uri) %>% 
            str_remove('.txt$') %>% 
            str_c('.rds') %>% 
            file.path(cache_dir, .)
        
        hgnc_complete <- 
            suppressWarnings(
                (function() read_tsv(hgnc_complete_uri, col_types = cols())) %>% 
                    cache(fn)) %>% 
            select(hgnc_id, symbol, name, location, ensembl_gene_id, entrez_id, alias_symbol, prev_symbol) %>% 
            mutate(entrez_id = as.integer(entrez_id))
        
        cavalier_cache$hgnc_complete <- hgnc_complete
    }
    return(hgnc_complete)
}

#' @importFrom tidyr replace_na separate_rows
#' @importFrom stringr str_c
#' @importFrom dplyr "%>%" mutate rename select if_else add_count filter case_when
get_hgnc_alias <- function() 
{
    hgnc_alias <- cavalier_cache$hgnc_alias
    
    if (is.null(hgnc_alias)) {
        
        hgnc_alias <-
            get_hgnc_complete() %>% 
            filter(!(is.na(alias_symbol) & is.na(prev_symbol))) %>%
            mutate(alias = case_when(
                !is.na(alias_symbol) & !is.na(prev_symbol) ~ str_c(alias_symbol, '|', prev_symbol),
                !is.na(alias_symbol)                       ~ alias_symbol,
                !is.na(prev_symbol)                        ~ prev_symbol)) %>% 
            select(symbol,  alias) %>% 
            separate_rows(alias, sep = '\\|') %>% 
            filter(! alias %in% symbol) %>% # remove ambiguities
            add_count(alias) %>%
            filter(n == 1) %>% # remove ambiguities
            select(-n) %>% 
            distinct() %>% 
            na.omit() %>% 
            arrange_all()
        
        cavalier_cache$hgnc_alias <- hgnc_alias
    }
    
    return(hgnc_alias)
}

#' @importFrom dplyr "%>%" select
get_hgnc_id <- function() 
{
    hgnc_id <- 
    
    if (is.null(hgnc_id)) {
        
        hgnc_id <-
            get_hgnc_complete() %>% 
            select(symbol, hgnc_id) %>% 
            na.omit()
        
        cavalier_cache$hgnc_id <- hgnc_id
    }
    
    return(hgnc_id)
}

#' @importFrom dplyr "%>%" select
get_hgnc_ensembl <- function() 
{
    hgnc_ensembl <- cavalier_cache$hgnc_ensembl
    
    if (is.null(hgnc_ensembl)) {
        
        hgnc_ensembl <-
            get_hgnc_complete() %>% 
            select(symbol, ensembl_gene_id) %>% 
            na.omit()
        
        cavalier_cache$hgnc_ensembl <- hgnc_ensembl
    }
    
    return(hgnc_ensembl)
}

#' @importFrom dplyr "%>%" select
get_hgnc_entrez <- function() 
{
    hgnc_entrez <- cavalier_cache$hgnc_entrez
    
    if (is.null(hgnc_entrez)) {
        
        hgnc_entrez <-
            get_hgnc_complete() %>% 
            select(symbol, entrez_id) %>% 
            na.omit()
        
        cavalier_cache$hgnc_entrez <- hgnc_entrez
    }
    
    return(hgnc_entrez)
}

# replace gene symbols with HGNC approved symbol
#' @export
hgnc_sym2sym <- function(symbols, remove_unknown = FALSE) {
    hgnc_alias <- get_hgnc_alias()
    at <- which(symbols %in% hgnc_alias$alias)
    ret <- replace(symbols, at, hgnc_alias$symbol[match(symbols[at], hgnc_alias$alias)])
    if (remove_unknown) { 
        ret[! ret %in% { get_hgnc_complete() %>% dplyr::pull(symbol) } ] <- NA_character_
    }
    return(ret)
}

#' @export
hgnc_sym2ensembl <- function(symbols) {
    symbols <- hgnc_sym2sym(symbols)
    hgnc_ensembl <- get_hgnc_ensembl()
    hgnc_ensembl$ensembl_gene_id[match(symbols, hgnc_ensembl$symbol)]
}

#' @export
hgnc_ensembl2sym <- function(ensembl_gene_id) {
    hgnc_ensembl <- get_hgnc_ensembl()
    hgnc_ensembl$symbol[match(ensembl_gene_id, hgnc_ensembl$ensembl_gene_id)]
}

#' @export
hgnc_id2sym <- function(ids) {
    hgnc_id <- get_hgnc_id()
    hgnc_id$symbol[match(ids, hgnc_id$hgnc_id)]
}

#' @export
hgnc_sym2entrez <- function(symbols) {
    symbols <- hgnc_sym2sym(symbols)
    hgnc_entrez <- get_hgnc_entrez()
    hgnc_entrez$entrez_id[match(symbols, hgnc_entrez$symbol)]
}

#' @export
hgnc_entrez2sym <- function(entrez_id) {
    hgnc_entrez <- get_hgnc_entrez()
    hgnc_entrez$symbol[match(entrez_id, hgnc_entrez$entrez_id)]
}
