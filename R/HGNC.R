
#' Get HGNC version from either get_cavalier_opt("hgnc_monthly_base_url") or disk cache
#' @importFrom dplyr slice pull mutate filter
get_hgnc_version <- function(
    db_mode = get_cavalier_opt("database_mode"),
    ver = get_cavalier_opt("hgnc_ver"))
{
  if (!is.null(ver)) {
    return(ver)
  }
  
  func_online <- function() {
    rvest::read_html('https://storage.googleapis.com/public-download-files') %>% 
      rvest::html_elements("key") %>%
      rvest::html_text() %>% 
      keep(str_detect, 'monthly/tsv/hgnc_complete_set_.+\\.txt$') %>% 
      str_extract('(?<=hgnc_complete_set_)\\d{4}-\\d{2}-\\d{2}') %>% 
      sort() %>% 
      last()
        
  }
  
  get_version(
    resource_name  = 'HGNC',
    cache_name = 'hgnc_complete_set',
    cache_subdir = 'HGNC',
    func_online = func_online,
    db_mode = db_mode
  )
}

#' Get HGNC complete table from either get_cavalier_opt("hgnc_monthly_base_url") or disk cache
#' @importFrom readr cols read_tsv
#' @importFrom dplyr "%>%" mutate rename
#' @importFrom stringr str_remove
#' @importFrom rlang is_scalar_character
get_hgnc_complete <- function(
    local_file = get_cavalier_opt("hgnc_local_file"),
    ver = get_hgnc_version()
)
{
    assert_that(
      is_scalar_character(ver),
      is_scalar_character(local_file) | is.null(local_file)
    )
    
    
    if (ver == "local") {
      con <- local_file
    } else {
      con <- str_c(get_cavalier_opt("hgnc_monthly_base_url"), 'hgnc_complete_set_', ver,'.txt')
    }
    
    fun <- function() {
      suppressWarnings(
        read_tsv(con, 
                 col_types = cols()
        )) %>% 
        select(hgnc_id, symbol, name, locus_group, location, ensembl_gene_id, entrez_id, alias_symbol, prev_symbol) %>% 
        mutate(entrez_id = as.integer(entrez_id))
    }
    
    if (ver == 'local') {
      return(
        fun()
      )
    } else {
      return(
        cache(
          fun = fun,
          name = 'hgnc_complete_set',
          version = ver,
          subdir = 'HGNC'
        )
      )
    }
}

#' Get HGNC gene symbols alias table to convert symbols to current HGNC symbol
#' @importFrom tidyr replace_na separate_rows
#' @importFrom stringr str_c
#' @importFrom dplyr arrange_all "%>%" mutate rename select if_else add_count filter case_when
get_hgnc_alias <- function() 
{
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
    arrange_all()
}

#' Get HGNC hgnc_id, symbol table
#' @importFrom dplyr "%>%" select distinct
get_hgnc_symbol <- function() 
{
  get_hgnc_complete() %>% 
    select(hgnc_id, symbol) %>% 
    distinct() %>% 
    na.omit()
}

#' Get HGNC hgnc_id, ensemble_gene_id table
#' @importFrom dplyr "%>%" select distinct
get_hgnc_ensembl <- function() 
{
    get_hgnc_complete() %>% 
     select(hgnc_id, ensembl_gene_id) %>% 
     distinct() %>% 
     na.omit()
}

#' Get HGNC hgnc_id, entrez_id table
#' @importFrom dplyr "%>%" select distinct
get_hgnc_entrez <- function() 
{
    get_hgnc_complete() %>% 
     select(hgnc_id, entrez_id) %>% 
     distinct() %>% 
     na.omit()
}

#' Get HGNC hgnc_id, locus_group table
#' @importFrom dplyr "%>%" select distinct
get_hgnc_locus_group <- function() 
{
  get_hgnc_complete() %>% 
    select(hgnc_id, locus_group) %>% 
    distinct() %>% 
    na.omit()
}

# Replace gene symbols with HGNC approved symbol
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

#' Convert hgnc_id to symbol
#' @export
hgnc_id2sym <- function(ids) 
{
  hgnc_symbol <- get_hgnc_symbol()
  hgnc_symbol$symbol[match(ids, hgnc_symbol$hgnc_id)]
}

#' Convert gene symbol to hgnc_id
#' @export
hgnc_sym2id <- function(symbols) 
{
  symbols <- hgnc_sym2sym(symbols)
  hgnc_symbol <- get_hgnc_symbol()
  hgnc_symbol$hgnc_id[match(symbols, hgnc_symbol$symbol)]
}

#' Convert hgnc_id to ensembl_gene_id
#' @export
hgnc_id2ensembl<- function(hgnc_ids) 
{
  hgnc_ensembl <- get_hgnc_ensembl()
  hgnc_ensembl$ensembl_gene_id[match(hgnc_ids, hgnc_ensembl$hgnc_id)]
}

#' Convert ensembl_gene_id to hgnc_id
#' @export
hgnc_ensembl2id <- function(ensembl_gene_ids) 
{
  hgnc_ensembl <- get_hgnc_ensembl()
  hgnc_ensembl$hgnc_id[match(ensembl_gene_ids, hgnc_ensembl$ensembl_gene_id)]
}

#' Convert hgnc_id to entrez_id
#' @export
hgnc_id2entrez<- function(hgnc_ids) 
{
  hgnc_entrez <- get_hgnc_entrez()
  hgnc_entrez$entrez_id[match(hgnc_ids, hgnc_entrez$hgnc_id)]
}

#' Convert entrez_id to hgnc_id
#' @export
hgnc_entrez2id <- function(entrez_ids) 
{
  hgnc_entrez <- get_hgnc_entrez()
  hgnc_entrez$hgnc_id[match(entrez_ids, hgnc_entrez$entrez_id)]
}

#' Convert entrez_id to symbol
#' @export
hgnc_entrez2sym <- function(entrez_ids) 
{
  hgnc_id2sym(hgnc_entrez2id(entrez_ids))
}

#' Convert symbol to entrez_id
#' @export
hgnc_sym2entrez <- function(symbols) 
{
  hgnc_id2entrez(hgnc_sym2id(symbols))
}

#' Convert ensembl_gene_id to symbol
#' @export
hgnc_ensembl2sym <- function(ensembl_gene_ids) 
{
  hgnc_id2sym(hgnc_ensembl2id(ensembl_gene_ids))
}

#' Convert symbol to ensembl_gene_id
#' @export
hgnc_sym2ensembl <- function(symbols) 
{
  hgnc_id2ensembl(hgnc_sym2id(symbols))
}

#' Convert ensembl_gene_id to entrez_id
#' @export
hgnc_ensembl2entrez <- function(ensembl_gene_ids) 
{
  hgnc_id2entrez(hgnc_ensembl2id(ensembl_gene_ids))
}

#' Convert entrez_id to ensembl_gene_id
#' @export
hgnc_entrez2ensembl <- function(entrez_ids) 
{
  hgnc_id2ensembl(hgnc_entrez2id(entrez_ids))
}

#' Get list of from HGNC by locus_group
get_hgnc_locus_group_list <- function(
    locus_group = c('protein-coding gene', 'non-coding RNA', 'pseudogene', 'other', 'ALL')
) 
{
  locus_group <- match.arg(locus_group)
  
  get_hgnc_locus_group() %>% 
    filter(!!locus_group == 'ALL' | locus_group == !!locus_group) %>% 
    select(-locus_group) %>% 
    left_join(get_hgnc_symbol(), by = 'hgnc_id') %>% 
    left_join(get_hgnc_ensembl(), by = 'hgnc_id') %>% 
    left_join(get_hgnc_entrez(), by = 'hgnc_id') %>% 
    mutate(list_id = str_c('HGNC:', locus_group),
           list_name = list_id,
           list_version = get_hgnc_version(),
           .before = 1)
}