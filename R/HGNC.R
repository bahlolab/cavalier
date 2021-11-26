
hgnc_complete_base_url <- 'http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/monthly/tsv/'

#' @importFrom dplyr slice pull mutate filter
get_hgnc_latest_version <- function()
{
    (function()
        tibble(filename = XML::getHTMLLinks(hgnc_complete_base_url)) %>% 
         filter(str_starts(filename, 'hgnc_complete_set_')) %>% 
         mutate(date = str_extract(filename, '(?<=hgnc_complete_set_)\\d{4}-\\d{2}-\\d{2}') %>% 
                    lubridate::as_date()) %>% 
         slice(which.max(date)) %>% 
         pull(date) %>% 
         as.character()) %>% 
        cache('hgnc_latest_version')
}

#' @importFrom readr cols read_tsv
#' @importFrom dplyr "%>%" mutate rename
#' @importFrom stringr str_remove
#' @importFrom rlang is_scalar_character
get_hgnc_complete <- function()
{
    ver <- get_hgnc_latest_version()
    uri <- str_c(hgnc_complete_base_url, 'hgnc_complete_set_', ver,'.txt')
    
    (function()
        suppressWarnings(read_tsv(uri, col_types = cols())) %>% 
            select(hgnc_id, symbol, name, location, ensembl_gene_id, entrez_id, alias_symbol, prev_symbol) %>% 
            mutate(entrez_id = as.integer(entrez_id))) %>% 
        cache(str_remove(basename(uri), '.txt$'),
              disk = TRUE)
}

#' @importFrom tidyr replace_na separate_rows
#' @importFrom stringr str_c
#' @importFrom dplyr arrange_all "%>%" mutate rename select if_else add_count filter case_when
get_hgnc_alias <- function() 
{
    (function()
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
         arrange_all()) %>% 
        cache('hgnc_alias')
}

#' @importFrom dplyr "%>%" select
get_hgnc_id <- function() 
{
    (function()
        get_hgnc_complete() %>% 
         select(symbol, hgnc_id) %>% 
         na.omit()) %>% 
        cache('hgnc_id')
}

#' @importFrom dplyr "%>%" select
get_hgnc_ensembl <- function() 
{
    (function()
        get_hgnc_complete() %>% 
         select(symbol, ensembl_gene_id) %>% 
         na.omit()) %>% 
        cache('hgnc_ensembl')
}

#' @importFrom dplyr "%>%" select
get_hgnc_entrez <- function() 
{
    (function()
        get_hgnc_complete() %>% 
         select(symbol, entrez_id) %>% 
         na.omit()) %>% 
        cache('hgnc_entrez')
}

# replace gene symbols with HGNC approved symbol
#' @export
hgnc_sym2sym <- function(symbols, remove_unknown = FALSE) 
{
    hgnc_alias <- get_hgnc_alias()
    at <- which(symbols %in% hgnc_alias$alias)
    ret <- replace(symbols, at, hgnc_alias$symbol[match(symbols[at], hgnc_alias$alias)])
    if (remove_unknown) { 
        ret[! ret %in% { get_hgnc_complete() %>% dplyr::pull(symbol) } ] <- NA_character_
    }
    return(ret)
}

#' @export
hgnc_sym2ensembl <- function(symbols) 
{
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
hgnc_id2sym <- function(ids) 
{
    hgnc_id <- get_hgnc_id()
    hgnc_id$symbol[match(ids, hgnc_id$hgnc_id)]
}

#' @export
hgnc_sym2entrez <- function(symbols) 
{
    symbols <- hgnc_sym2sym(symbols)
    hgnc_entrez <- get_hgnc_entrez()
    hgnc_entrez$entrez_id[match(symbols, hgnc_entrez$symbol)]
}

#' @export
hgnc_entrez2sym <- function(entrez_id) 
{
    hgnc_entrez <- get_hgnc_entrez()
    hgnc_entrez$symbol[match(entrez_id, hgnc_entrez$entrez_id)]
}
