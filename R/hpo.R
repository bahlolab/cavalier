gene_2_phen <- function() 
{
  col_names <- 
    c('entrez_gene_id', 'entrez_gene_symbol', 'hpo_term_id', 'hpo_term_name', 'frequency_raw', 
      'frequency_hpo', 'additional_info_from_g_d_source', 'g_d_source', 'disease_id_for_link')
  url <- 'http://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt'  
  url <- 'https://ci.monarchinitiative.org/view/hpo/job/hpo.annotations/lastSuccessfulBuild/artifact/rare-diseases/util/annotation/genes_to_phenotype.txt'
  
  g2p <-
    readr::read_tsv(url, col_names = col_names, skip = 1)  %>% 
    mutate(entrez_gene_id = as.integer(entrez_gene_id)) %>% 
    select(entrez_gene_id, entrez_gene_symbol, hpo_term_id, hpo_term_name, disease_id=disease_id_for_link) %>% 
    (function(x)
      left_join(
        anti_join(x, 
                  get_inheritance_terms() %>% mutate(inheritance = coalesce(group, name)),
                  by = c(hpo_term_id = 'ontologyId')), 
        #inheritance
        inner_join(x,
                   get_inheritance_terms() %>% mutate(inheritance = coalesce(group, name)),
                   by = c(hpo_term_id = 'ontologyId')) %>% 
          group_by(entrez_gene_id, entrez_gene_symbol, disease_id) %>% 
          summarise(inheritance = str_c(sort(inheritance), collapse = ' & '),
                    .groups = 'drop'),
        by = c("entrez_gene_id", "entrez_gene_symbol", "disease_id"))
     )
}

get_hpo_gene_list2 <- function(hpo_id)
{
  g2p %>% 
    filter(hpo_term_id == hpo_id) %>% 
    select(entrez_gene_id, disease_id_for_link)
}