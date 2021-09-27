hpo_jenkins_url <- 'https://ci.monarchinitiative.org/view/hpo/job/hpo.annotations/lastSuccessfulBuild'
hpo_jenkins_g2p_url <- 'https://ci.monarchinitiative.org/view/hpo/job/hpo.annotations/lastSuccessfulBuild/artifact/rare-diseases/util/annotation/genes_to_phenotype.txt'
hpo_jenkins_p2g_url <- 'https://ci.monarchinitiative.org/view/hpo/job/hpo.annotations/lastSuccessfulBuild/artifact/rare-diseases/util/annotation/phenotype_to_genes.txt'


get_latest_hpo_build_num <- function()
{
  content(GET(str_c(hpo_jenkins_url, '/buildNumber')))
}

# this contains childmost term for a phenotype heirachy
# use for annotating genes
get_genes_to_phenotype <- function() 
{
  genes_to_phenotype <- cavalier_cache$genes_to_phenotype
  
  if (is.null(genes_to_phenotype)) {
    
    col_names <- 
      c('entrez_gene_id', 'entrez_gene_symbol', 'hpo_term_id', 'hpo_term_name', 'frequency_raw', 
        'frequency_hpo', 'additional_info', 'g_d_source', 'disease_id')
    
    build_num <- get_latest_hpo_build_num()
    
    genes_to_phenotype <- 
      (function() read_tsv(hpo_jenkins_g2p_url,
                           col_names = col_names,
                           skip = 1,
                           col_types = cols())) %>% 
      cache(str_replace(basename(hpo_jenkins_g2p_url),
                        '\\.txt', str_c('_v', build_num, '.txt')))
    
    cavalier_cache$genes_to_phenotype <- genes_to_phenotype
  }
  
  return(genes_to_phenotype)  
}

# this contains all hpo terms
# use for creating gene lists from hpo terms
get_phenotype_to_genes<- function() 
{
  phenotype_to_genes <- cavalier_cache$phenotype_to_genes
  
  if (is.null(phenotype_to_genes)) {
    
    col_names <- 
      c('hpo_term_id', 'hpo_term_name', 'entrez_gene_id', 'entrez_gene_symbol',
        'additional_info', 'g_d_source', 'disease_id')
    
    build_num <- get_latest_hpo_build_num()
    
    phenotype_to_genes <- 
      (function() read_tsv(hpo_jenkins_p2g_url,
                           col_names = col_names,
                           skip = 1,
                           col_types = cols())) %>% 
      cache(str_replace(basename(hpo_jenkins_p2g_url),
                        '\\.txt', str_c('_v', build_num, '.txt')))
    
    cavalier_cache$phenotype_to_genes <- phenotype_to_genes
  }
  
  return(phenotype_to_genes)  
}

get_term_names <- function() 
{
  term_names <- cavalier_cache$term_names
  
  if (is.null(term_names)) {
    
    term_names <-
      get_phenotype_to_genes() %>% 
      select(hpo_term_id, hpo_term_name) %>% 
      distinct()
    
    cavalier_cache$term_names <- term_names
  }
  
  return(term_names)  
}

gene_2_phen <- function() 
{
  
  
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
