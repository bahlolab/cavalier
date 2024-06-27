
.onLoad <- function(libname, pkgname) {
  
  # create memory-cached version of cavalier functions to save time
  
  use_memoisation <- get_cavalier_opt("use_memoisation")
  
  if (!is.null(use_memoisation) && use_memoisation) {
    
    ### HPO ###
    get_hpo_version <<- memoise::memoise(get_hpo_version)
    hpo_api_get <<- memoise::memoise(hpo_api_get)
    get_hpo_g2p_p2g <<- memoise::memoise(get_hpo_g2p_p2g)
    get_gene_disease_map <<- memoise::memoise(get_gene_disease_map)
    get_hpo_term_names <<- memoise::memoise(get_hpo_term_names)
    get_hpo_gene_list <<- memoise::memoise(get_hpo_gene_list)
    get_disease_name_cache <<- memoise::memoise(get_disease_name_cache)
    
    ### HGNC ###
    get_hgnc_version <<- memoise::memoise(get_hgnc_version)
    get_hgnc_complete <<- memoise::memoise(get_hgnc_complete)
    get_hgnc_alias <<- memoise::memoise(get_hgnc_alias)
    get_hgnc_symbol <<- memoise::memoise(get_hgnc_symbol)
    get_hgnc_ensembl <<- memoise::memoise(get_hgnc_ensembl)
    get_hgnc_entrez <<- memoise::memoise(get_hgnc_entrez)
    get_hgnc_locus_group <<- memoise::memoise(get_hgnc_locus_group)
    get_hgnc_locus_group_list <<- memoise::memoise(get_hgnc_locus_group_list)
    
    ### Gene intolerance ###
    get_gevir_table <<- memoise::memoise(get_gevir_table)
    get_gevir_table_hgnc <<- memoise::memoise(get_gevir_table_hgnc)
    
    ### PanelApp ###
    get_panelapp_panels <<- memoise::memoise(get_panelapp_panels)
    get_panelapp_panel_version <<- memoise::memoise(get_panelapp_panel_version)
    get_panelapp_panel <<- memoise::memoise(get_panelapp_panel)
    
    ### GTEx ###
    get_gtex_expression_hgnc <<- memoise::memoise(get_gtex_expression_hgnc)
    
    ### UCSC ###
    get_centromeres_gaps <<- memoise::memoise(get_centromeres_gaps)
    
    ### Genes4Epilepsy
    get_g4e_version <<- memoise::memoise(get_g4e_version)
    get_g4e_full_list <<- memoise::memoise(get_g4e_full_list)
  }
  
  # set cavalier options from global options
  options() %>% 
    (function(x) {
      x <- x[str_starts(names(x), 'cavalier.')]
      names(x) <- str_remove(names(x), '^cavalier.')
      x
    }) %>% 
    (function(x) do.call(set_cavalier_opt, x))
    

  invisible()
  
}
