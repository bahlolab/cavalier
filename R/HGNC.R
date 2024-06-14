#### hgnc_ids used as primary identifier when joining different gene identifiers

#' @importFrom dplyr slice pull mutate filter
get_hgnc_latest_version <- function()
{
  hgnc_monthly_base_url <- get_cavalier_opt("hgnc_monthly_base_url")
  (function()
    tryCatch(
      retry('GET', hgnc_monthly_base_url) %>% 
        content(encoding = 'UTF-8') %>% 
        rvest::html_nodes('a') %>%
        rvest::html_attr("href") %>% 
        tibble(filename = .) %>% 
        filter(str_starts(filename, 'hgnc_complete_set_')) %>% 
        mutate(date = str_extract(filename, '(?<=hgnc_complete_set_)\\d{4}-\\d{2}-\\d{2}') %>% 
                 lubridate::as_date()) %>% 
        slice(which.max(date)) %>% 
        pull(date) %>% 
        as.character(),
      error = function(e) {
        ver <-
          tibble(filename = list.files(get_cache_dir(), pattern = '^hgnc_complete_set_\\d{4}-\\d{2}-\\d{2}.*\\.rds$')) %>% 
          filter(str_starts(filename, 'hgnc_complete_set_')) %>% 
          mutate(date = str_extract(filename, '(?<=hgnc_complete_set_)\\d{4}-\\d{2}-\\d{2}') %>% 
                   lubridate::as_date()) %>% 
          slice(which.max(date)) %>% 
          pull(date) %>% 
          as.character()
        hgnc_files <- 
        if (length(ver)) {
          warning("Coudn't access latest HGNC build at: ", hgnc_monthly_base_url, '. ',
                  "Using cached version ", ver, '.')
          ver
        } else {
          stop('could not get HGNC latest version')
        }
      }
    )) %>% 
    cache('hgnc_latest_version')
}

#' @importFrom readr cols read_tsv
#' @importFrom dplyr "%>%" mutate rename
#' @importFrom stringr str_remove
#' @importFrom rlang is_scalar_character
get_hgnc_complete <- function()
{
    ver <- get_cavalier_opt("hgnc_ver")
    assert_that(is_scalar_character(ver))
    
    if (ver == "latest") {
      ver <- get_hgnc_latest_version()
    }
    if (ver == "local") {
      con <- get_cavalier_opt("hgnc_local_file")
    } else {
      con <- str_c(get_cavalier_opt("hgnc_monthly_base_url"), 'hgnc_complete_set_', ver,'.txt')
    }
    
    cache_name <- basename(con)
    
    fun <- function() {
      suppressWarnings(
        read_tsv(con, 
                 col_types = cols()
        )) %>% 
        select(hgnc_id, symbol, name, location, ensembl_gene_id, entrez_id, alias_symbol, prev_symbol) %>% 
        mutate(entrez_id = as.integer(entrez_id))
    }
    
    cache(
      fun = fun,
      name = cache_name,
      disk = ver != "local",
      subdir = 'HGNC')
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
         filter(! alias %in% get_hgnc_complete()$symbol) %>% # remove ambiguities
         add_count(alias) %>%
         filter(n == 1) %>% # remove ambiguities
         select(-n) %>% 
         distinct() %>% 
         na.omit() %>% 
         arrange_all()) %>% 
        cache('hgnc_alias')
}

#' @importFrom dplyr "%>%" select distinct
get_hgnc_symbol <- function() 
{
    (function()
        get_hgnc_complete() %>% 
         select(hgnc_id, symbol) %>% 
         distinct() %>% 
         na.omit()) %>% 
        cache('hgnc_symbol')
}

#' @importFrom dplyr "%>%" select distinct
get_hgnc_ensembl <- function() 
{
    (function()
        get_hgnc_complete() %>% 
         select(hgnc_id, ensembl_gene_id) %>% 
         distinct() %>% 
         na.omit()) %>% 
        cache('hgnc_ensembl')
}

#' @importFrom dplyr "%>%" select distinct
get_hgnc_entrez <- function() 
{
    (function()
        get_hgnc_complete() %>% 
         select(hgnc_id, entrez_id) %>% 
         distinct() %>% 
         na.omit()) %>% 
        cache('hgnc_entrez')
}

# replace gene symbols with HGNC approved symbol
#' @export
hgnc_sym2sym <- function(symbols, remove_unknown = FALSE) 
{
    unknown <- which(!symbols %in% get_hgnc_complete()$symbol)
    hgnc_alias <- get_hgnc_alias()
    at <- unknown[symbols[unknown] %in% hgnc_alias$alias]
    ret <- replace(symbols, at, hgnc_alias$symbol[match(symbols[at], hgnc_alias$alias)])
    if (remove_unknown) { 
        at <- setdiff(unknown, at)
        ret[at] <- NA_character_
    }
    return(ret)
}

#' @export
hgnc_id2sym <- function(ids) 
{
  hgnc_symbol <- get_hgnc_symbol()
  hgnc_symbol$symbol[match(ids, hgnc_symbol$hgnc_id)]
}

#' @export
hgnc_sym2id <- function(symbols) 
{
  symbols <- hgnc_sym2sym(symbols)
  hgnc_symbol <- get_hgnc_symbol()
  hgnc_symbol$hgnc_id[match(symbols, hgnc_symbol$symbol)]
}

#' @export
hgnc_id2ensembl<- function(hgnc_ids) 
{
  hgnc_ensembl <- get_hgnc_ensembl()
  hgnc_ensembl$ensembl_gene_id[match(hgnc_ids, hgnc_ensembl$hgnc_id)]
}

#' @export
hgnc_ensembl2id <- function(ensembl_gene_ids) 
{
  hgnc_ensembl <- get_hgnc_ensembl()
  hgnc_ensembl$hgnc_id[match(ensembl_gene_ids, hgnc_ensembl$ensembl_gene_id)]
}

#' @export
hgnc_id2entrez<- function(hgnc_ids) 
{
  hgnc_entrez <- get_hgnc_entrez()
  hgnc_entrez$entrez_id[match(hgnc_ids, hgnc_entrez$hgnc_id)]
}

#' @export
hgnc_entrez2id <- function(entrez_ids) 
{
  hgnc_entrez <- get_hgnc_entrez()
  hgnc_entrez$hgnc_id[match(entrez_ids, hgnc_entrez$entrez_id)]
}

#' @export
hgnc_entrez2sym <- function(entrez_ids) 
{
  hgnc_id2sym(hgnc_entrez2id(entrez_ids))
}

#' @export
hgnc_sym2entrez <- function(symbols) 
{
  hgnc_id2entrez(hgnc_sym2id(symbols))
}

#' @export
hgnc_ensembl2sym <- function(ensembl_gene_ids) 
{
  hgnc_id2sym(hgnc_ensembl2id(ensembl_gene_ids))
}

#' @export
hgnc_sym2ensembl <- function(symbols) 
{
  hgnc_id2ensembl(hgnc_sym2id(symbols))
}
