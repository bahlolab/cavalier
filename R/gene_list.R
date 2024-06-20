
#' Return a gene list from PanelApp, HPO or Genes4Epilepsy
#' 
#' Return data frame with columns list_id, list_name, list_version, symbol, hgnc_id, ensembl_gene_id, entrez_id, inheritance
#' NB confidence_level only used for PanelApp.
#' @examples
#' # Get PanelApp Australia Panel 202 (Genetic Epilespsy)
#' get_gene_list('PAA:202', version = '1.26')
#' # Get HPO genes associated with HP:0001250 (Seizure)
#' get_gene_list('HP:0001250')
#' #' Get all Genes4Epilepsy genes
#' get_gene_list('G4E:Epilepsy')
#' # Get Genes4Epilepsy genes associated Malformations of Cortical Development
#' get_gene_list('G4E:MCD')
#' @export
get_gene_list <- function(id, version = NULL, save = NULL, min_confidence = 2L)
{
  assert_that(is_scalar_character(id))
  
  ## assume PanelApp starts with PA
  pref <- str_extract(id, '^[^:]+')
  
  if (pref == 'HP') {
    result <- 
      get_hpo_gene_list(id, hpo_version = version) %>% 
      mutate(hgnc_id = hgnc_entrez2id(entrez_id),
             ensembl_gene_id = hgnc_id2ensembl(hgnc_id))
  } else if (str_starts(pref, 'PA')) {
    result <-
      get_panelapp_panel(id, version = version) %>%
      filter(confidence_level > min_confidence) %>% 
      mutate(symbol = coalesce(hgnc_id2sym(hgnc_id), hgnc_sym2sym(gene)),
             ensembl_gene_id = hgnc_id2ensembl(hgnc_id),
             entrez_id = hgnc_id2entrez(hgnc_id)) %>% 
      rename(list_version = version)
      
  } else if (pref == 'G4E') {
    pheno <- str_extract(id, '(?<=:).+')
    result <-
      get_g4e_phenotype_list(version = version, phenotype = pheno) %>% 
      mutate(list_name = str_c('Genes4Epilepsy - ', pheno),
             symbol = coalesce(hgnc_id2sym(hgnc_id), hgnc_sym2sym(symbol)))
  } else {
    stop('Unrecognised web gene list id "', id, '"')
  }
  
  result <-
    result %>% 
    select(
      list_id,
      list_name,
      list_version,
      symbol, 
      hgnc_id,
      ensembl_gene_id,
      entrez_id,
      inheritance
    ) %>% 
    arrange_all() %>% 
    distinct()
  
  if (!is.null(save)) {
    write_tsv(result, save)
    invisible(result)
  } else {
    result
  }
}


#' @export
get_gene_list_versions <- function(ids, save = NULL)
{
  assert_that(is_character(ids))
  
  result <-
    tibble(id = ids) %>% 
    distinct() %>% 
    mutate(version = map_chr(id, function(id) {
      get_gene_list(id) %>% 
        pull(list_version) %>% 
        first()
    }))
  
  if (!is.null(save)) {
    write_tsv(result, save)
    invisible(result)
  } else {
    result
  }
}
