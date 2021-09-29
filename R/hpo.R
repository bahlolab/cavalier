hpo_jenkins_url <- 'https://ci.monarchinitiative.org/view/hpo/job/hpo.annotations/lastSuccessfulBuild'
hpo_jenkins_g2p_url <- 'https://ci.monarchinitiative.org/view/hpo/job/hpo.annotations/lastSuccessfulBuild/artifact/rare-diseases/util/annotation/genes_to_phenotype.txt'
hpo_jenkins_p2g_url <- 'https://ci.monarchinitiative.org/view/hpo/job/hpo.annotations/lastSuccessfulBuild/artifact/rare-diseases/util/annotation/phenotype_to_genes.txt'
inheritance_term_id <- 'HP:0000005'

## TODO: move these to options in exported functions(e.g. secure = TRUE), use as httr::with_config
#'@export
insecure <- function() httr::set_config(httr::config(ssl_verifypeer = 0L))
#'@export
secure <- function() httr::set_config(httr::config(ssl_verifypeer = 1L))


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
get_phenotype_to_genes <- function() 
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

#'@importFrom httr GET accept_json content
#'@importFrom stringr str_c
hpo_api_get <- function(extension,
                        base_url = "https://hpo.jax.org/api/hpo/")
{
  url <- str_c(base_url, extension)
  
  response <- GET(url, accept_json())
  
  if (response$status_code == 200) {
    return(content(response))
  }
  return(NULL)
}

#'@importFrom dplyr tibble as_tibble bind_rows "%>%"
#'@importFrom urltools url_encode
#'@importFrom memoise memoise 
#'@importFrom purrr map_df
#'@export
get_hpo_term <- function(hpo_id)
{
  assert_that(is.character(hpo_id),
              all(str_detect(hpo_id, 'HP:\\d{7}'), na.rm = TRUE))
  
  mapper <- cavalier_cache$hpo_get_term_mapper
  
  if (is.null(mapper)) {
    mapper <- memoise(
      function(hpo_id) {
        if (!is.na(hpo_id)) {
          result <- hpo_api_get(str_c('term/', url_encode(hpo_id)))
          if (!is.null(result)) {
            data <-
              result$details %>%
              map(~ `if`(is.list(.), list(unlist(.)), .)) %>% 
              as_tibble() %>% 
              mutate(parents = list(bind_rows(result$relations$parents)),
                     children = list(bind_rows(result$relations$children))) %>% 
              rename(ontologyId = id) %>% 
              select(ontologyId, everything())
            return(data)
          } 
        }
        return(tibble(ontologyId = hpo_id))
      })
    cavalier_cache$hpo_get_term_mapper <- mapper
  }
  
  map_df(hpo_id, mapper)
}

#'@importFrom rlang is_integerish 
#'@importFrom urltools url_encode
#'@export
get_hpo_term_genes <- function(hpo_id)
{
  assert_that(is.character(hpo_id),
              all(str_detect(hpo_id, 'HP:\\d{7}'), na.rm = TRUE))
  
  mapper <- cavalier_cache$get_hpo_term_genes_mapper
  
  if (is.null(mapper)) {
    mapper <- memoise(
      function(hpo_id) {
        if (!is.na(hpo_id)) {
          result <- hpo_api_get(str_c('term/', url_encode(hpo_id[[1]]), '/genes?max=-1'))
          if (!is.null(result)) {
            data <-
              tibble(ontologyId = hpo_id,
                     genes = list(
                       map_df(result$genes,
                              ~ tibble(entrezGeneId = .$entrezGeneId,
                                       entrezGeneSymbol = .$entrezGeneSymbol,
                                       dbDiseases = list(bind_rows(.$dbDiseases))))
                     ))
            return(data)
          } 
        }
        return(tibble(ontologyId = hpo_id))
      })
    cavalier_cache$get_hpo_term_genes_mapper <- mapper
  }
  
  map_df(hpo_id, mapper)
}

#'@importFrom urltools url_encode
#'@importFrom rlang is_integerish 
#'@export
get_hpo_gene <- function(entrez_id)
{
  assert_that(is_integerish(entrez_id),
              all(entrez_id > 0, na.rm = TRUE))
  
  mapper <- cavalier_cache$get_hpo_gene_mapper
  
  if (is.null(mapper)) {
    mapper <- memoise(
      function(entrez_id) {
        if (!is.na(entrez_id)) {
          result <- hpo_api_get(str_c('gene/', entrez_id))
          if (!is.null(result)) {
            data <- 
              as_tibble(result$gene) %>% 
              mutate(termAssoc = list(bind_rows(result$termAssoc)),
                     diseaseAssoc = list(bind_rows(result$diseaseAssoc)))
            return(data)
          } 
        }
        return(tibble(entrezGeneId = entrez_id))
      })
    cavalier_cache$get_hpo_gene_mapper <- mapper
  }
  
  # map_df_prog(entrez_id, mapper)
  map_df(entrez_id, mapper)
}

inheritance_groups <- c(
  'Autosomal dominant' = 'HP:0000006',
  'Autosomal recessive' = 'HP:0000007',
  'X-linked' = 'HP:0001417',
  'Y-linked' = 'HP:0001450')

get_inheritance_terms <- function() 
{
  
  inheritance_terms <- cavalier_cache$inheritance_terms
  terms_exclude <- 'HP:0001425'
  
  if (is.null(inheritance_terms)) {
    groups <- inheritance_groups
    
    inheritance_terms <- 
      get_hpo_term('HP:0000005') %>% 
      select(children) %>% 
      unnest(children) %>% 
      mutate(group = names(groups)[match(ontologyId, groups)]) %>% 
      (function(x) 
        filter(x, childrenCount > 0) %>% 
         select(group, ontologyId) %>%
         mutate(children = map(ontologyId, get_hpo_term)) %>%
         select(group, children) %>%
         unnest(children) %>%
         select(group, children) %>%
         unnest(children) %>% 
         mutate(group = coalesce(names(groups)[match(ontologyId, groups)],
                                 group)) %>% 
         (function(y)
           filter(y, childrenCount > 0) %>%
            select(group, ontologyId) %>%
            mutate(children = map(ontologyId, get_hpo_term)) %>%
            select(group, children) %>%
            unnest(children) %>%
            select(group, children) %>%
            unnest(children) %>% 
            bind_rows(y, .)
         ) %>% 
         bind_rows(x, .)
      ) %>% 
      select(ontologyId, name, group) %>% 
      # filter(!ontologyId %in% terms_exclude) %>% 
      arrange(ontologyId) %>% 
      mutate(name = str_remove(name, '\\sinheritance$'))
    
    cavalier_cache$inheritance_terms <- inheritance_terms
  }
  
  return(inheritance_terms)
  
}

#'@importFrom purrr map_int
get_term_descendants <- function(term_id) 
{
  recursion <- cavalier_cache$get_term_descendants_recursion
  
  if (is.null(recursion)) {
    recursion <- memoise(function(init) {
      init %>%
        filter(childrenCount > 0) %>%
        pull(ontologyId) %>%
        map_df(function(id) {
          get_hpo_term(id) %>% 
            select(children) %>% 
            unnest(children) %>% 
            mutate(parentId = !!id) %>% 
            select(ontologyId, name, parentId, childrenCount) %>% 
            recursion()
        }) %>% 
        bind_rows(init, .)
    })
    cavalier_cache$get_term_descendants_recursion <- recursion
  }
  
  get_hpo_term(term_id) %>% 
    mutate(parentId = NA_character_,
           childrenCount = map_int(children, nrow)) %>% 
    select(ontologyId, name, parentId, childrenCount) %>% 
    recursion()
}

#'@importFrom urltools url_encode
#'@export
get_hpo_disease <- function(disease_id)
{
  assert_that(
    is.character(disease_id),
    all(str_detect(disease_id, '^(OMIM)|(ORPHA):\\d+$'), na.rm = TRUE))
  
  mapper <- cavalier_cache$get_hpo_disease_mapper
  
  if (is.null(mapper)) {
    mapper <- memoise(
      function(disease_id) {
        if (!is.na(disease_id)) {
          result <- hpo_api_get(str_c('disease/', url_encode(disease_id)))
          if (!is.null(result)) {
            data <-
              as_tibble(result$disease) %>% 
              mutate(
                catTermsMap = list(map_df(result$catTermsMap, function(data) {
                  tibble(catLabel = data$catLabel,
                         terms = list(bind_rows(data$terms))) %>% 
                    unnest(terms)
                })),
                geneAssoc = list(bind_rows(result$geneAssoc)))
            return(data)
          } 
        }
        return(tibble(diseaseId = disease_id))
      })
    cavalier_cache$get_hpo_disease_mapper <- mapper
  }
  
  # map_df_prog(disease_id, mapper)
  map_df(disease_id, mapper)
}

#' @export
#' @importFrom dplyr inner_join anti_join add_row
get_hpo_gene_list <- function(hpo_id) {
  
  assert_that(is_scalar_character(hpo_id),
              !is.na(hpo_id),
              str_detect(hpo_id, 'HP:\\d{7}'))
  
  term_name <- 
    with(get_term_names(),
         hpo_term_name[match(hpo_id, hpo_term_id)])
  
  get_phenotype_to_genes() %>% 
    filter(hpo_term_id == hpo_id) %>% 
    select(entrez_gene_id, entrez_gene_symbol, disease_id) %>%
    (function(gd) {
      gd %>% 
        inner_join(get_genes_to_phenotype() %>% 
                     select(hpo_term_id, entrez_gene_id, disease_id),
                   by = c("entrez_gene_id", "disease_id")) %>% 
        inner_join(get_term_descendants(inheritance_term_id) %>% 
                     select(hpo_term_id = ontologyId),
                   by = "hpo_term_id") %>% 
        rename(inheritance = hpo_term_id) %>% 
        chop(inheritance) %>% 
        complete(gd) %>% 
        mutate(inheritance = simplify_inheritance(inheritance))
    }) %>% 
    select(gene = entrez_gene_symbol, entrez_gene_id, disease_id, inheritance) %>% 
    arrange(gene, disease_id) %>% 
    mutate(., version = str_c('build_', get_latest_hpo_build_num())) %>% 
    mutate(list_id = hpo_id,
           list_name = term_name) %>% 
    select(list_id, list_name, everything())
}

hpo_simple_inheritnace <-
  tibble(name = character(), ontologyId = list()) %>% 
  add_row(name = 'AR', ontologyId = list('HP:0000007')) %>% 
  add_row(name = 'AD', ontologyId = list('HP:0000006')) %>% 
  add_row(name = 'AD', ontologyId = list('HP:0025352')) %>% 
  add_row(name = 'XL', ontologyId = list('HP:0001417')) %>% 
  add_row(name = 'XL', ontologyId = list('HP:0001419')) %>% 
  add_row(name = 'XL', ontologyId = list('HP:0001423')) %>% 
  add_row(name = 'XL', ontologyId = list(c('HP:0001417', 'HP:0001419'))) %>% 
  add_row(name = 'XL', ontologyId = list(c('HP:0001417', 'HP:0001423'))) %>% 
  add_row(name = 'XL', ontologyId = list(c('HP:0001419', 'HP:0001423'))) %>% 
  add_row(name = 'MT', ontologyId = list('HP:0001427')) %>% 
  add_row(name = 'MT', ontologyId = list(c('HP:0000007', 'HP:0001427'))) %>% 
  add_row(name = 'AR/AD', ontologyId = list(c('HP:0000006', 'HP:0000007')))

hpo_inheritnace_ignore <- c(
  'HP:0001425', # Heterogeneous 
  'HP:0003745', # Sporadic,
  'HP:0003743', # Genetic anticipation
  'HP:0001428', # Somatic mutation
  'HP:0001426', # Multifactorial inheritance
  'HP:0001452', # Autosomal dominant contiguous gene syndrome
  'HP:0001466', # Contiguous gene syndrome
  'HP:0001442', # Somatic mosaicism
  'HP:0010984', # Digenic inheritance
  'HP:0010982', # Polygenic inheritance
  'HP:0003744'  # Genetic anticipation with paternal anticipation bias
)

simplify_inheritance <- function(inh)
{
  tibble(ontologyId = as.list(inh)) %>% 
    mutate(ontologyId = map(ontologyId, ~ sort(setdiff(., hpo_inheritnace_ignore)))) %>%
    left_join(hpo_simple_inheritnace, by = "ontologyId") %>% 
    mutate(name = replace_na(name, 'OTHER')) %>%
    pull(name)
}
