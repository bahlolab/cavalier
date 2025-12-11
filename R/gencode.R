
get_gencode_version_latest <- function() {
  base_url <- 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/'
  
  rvest::read_html(base_url) %>% 
    rvest::html_elements("a") %>% 
    rvest::html_attr("href") %>% 
    # Filter for the specific Primary Assembly GTF pattern
    keep(str_detect, 'gencode\\.v\\d+\\.primary_assembly\\.annotation\\.gtf\\.gz$') %>% 
    # In case there are duplicates (rare in latest_release), take the first match
    first() %>% 
    str_extract('(?<=gencode\\.)v[0-9]+')
}

get_gencode_version <- function(db_mode = get_cavalier_opt("database_mode"),
                                ver = get_cavalier_opt("gencode_version")) {
  if (is.null(ver)) {
    ver <-
      get_version(
        resource_name  = 'Gencode',
        cache_name = 'gencode_coords',
        cache_subdir = 'Gencode',
        func_online = get_gencode_version_latest,
        db_mode = db_mode
      )
  }
  
  message("Using Gencode version ", ver)
  
  ver
}

#' Get HGNC complete table from either get_cavalier_opt("hgnc_monthly_base_url") or disk cache
#' @importFrom rlang is_scalar_character
#' @export
get_gencode_coords <- function(
    ver = get_gencode_version()
)
{
    assert_that(
      is_scalar_character(ver)
    )
    
    fun <- function() {
      base_url <- 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/'
      url <- str_c(base_url, 'gencode.', ver, '.primary_assembly.annotation.gtf.gz')
      cmd <- str_c(
        "curl -s ",
        url, 
        " | gunzip -c | grep -P '\tgene\t' | ",
        "awk -F'\t' '{ match($9, /gene_id \"([^\"]+)\"/, gid); match($9, /hgnc_id \"([^\"]+)\"/, hid); print $1 \"\t\" $4 \"\t\" $5 \"\t\" gid[1] \"\t\" hid[1] }'"
      )
      
      read_tsv(
        pipe(cmd), 
        col_names = c("chromosome", "start", "end", "ensembl_gene_id", "hgnc_id"),
        col_types = "ciicc"
        ) %>% 
        mutate(ensembl_gene_id = str_remove(ensembl_gene_id, '\\.[0-9]+$'))
    }

    return(
      cache(
        fun = fun,
        name = 'gencode_coords',
        version = ver,
        subdir = 'Gencode'
      )
    )
}
