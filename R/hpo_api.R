
#'@export
insecure <- function() httr::set_config(httr::config(ssl_verifypeer = 0L))
#'@export
secure <- function() httr::set_config(httr::config(ssl_verifypeer = 1L))

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
      arrange(ontologyId) %>% 
      mutate(name = str_remove(name, '\\sinheritance$'))
    
    cavalier_cache$inheritance_terms <- inheritance_terms
  }
  
  return(inheritance_terms)
  
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
#' @importFrom dplyr inner_join anti_join
get_hpo_gene_list <- function(hpo_id) {
  
  assert_that(is_scalar_character(hpo_id),
              !is.na(hpo_id),
              str_detect(hpo_id, 'HP:\\d{7}'))
  
  get_hpo_term_genes(hpo_id) %>% 
    unnest(genes) %>% 
    select(entrezGeneId, entrezGeneSymbol) %>% 
    (function(genes)
      genes %>% 
       pull(entrezGeneId) %>% 
       get_hpo_gene() %>% 
       select(entrezGeneId, entrezGeneSymbol, termAssoc) %>% 
       unnest(termAssoc) %>% 
       select(entrezGeneId, entrezGeneSymbol, ontologyId) %>% 
       inner_join(get_inheritance_terms(), 'ontologyId') %>%
       add_count(entrezGeneId) %>% 
       (function(x)
         filter(x, n > 1) %>% 
          pull(entrezGeneId) %>% 
          unique() %>% 
          get_hpo_gene() %>% 
          select(entrezGeneId,entrezGeneSymbol, diseaseAssoc) %>%
          unnest(diseaseAssoc) %>%
          select(entrezGeneId, entrezGeneSymbol, diseaseId) %>% 
          (function(y)
            pull(y, diseaseId) %>% 
             unique() %>% 
             get_hpo_disease() %>% 
             select(diseaseId, catTermsMap) %>%
             unnest(catTermsMap) %>%
             group_by(diseaseId) %>% 
             filter(any(ontologyId == hpo_id)) %>% 
             ungroup() %>% 
             select(diseaseId, ontologyId) %>% 
             inner_join(get_inheritance_terms(), c('ontologyId')) %>%
             inner_join(y, ., by = 'diseaseId') %>% 
             select(-diseaseId)
          ) %>% 
          (function(y) 
            bind_rows(anti_join(x, y, by = 'entrezGeneId'), y)
          ) %>% 
          select(-n)
       ) %>% 
       complete(genes)
    ) %>% 
    arrange(entrezGeneId, ontologyId) %>% 
    mutate(inheritance = coalesce(group, name)) %>% 
    select(entrezGeneId, entrezGeneSymbol, inheritance) %>% 
    distinct() %>% 
    group_by(entrezGeneId, entrezGeneSymbol) %>% 
    summarise(inheritance =  str_c(inheritance, collapse = ' & '),
              .groups = 'drop') %>% 
    select(gene = entrezGeneSymbol,
           inheritance = inheritance) %>% 
    arrange(gene) %>% 
    mutate(., version = digest_df(.)) %>% 
    mutate(list_id = hpo_id,
           list_name = get_hpo_term(hpo_id)$name) %>% 
    select(list_id, list_name, everything())
}

