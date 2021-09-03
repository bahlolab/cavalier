
httr::set_config(httr::config(ssl_verifypeer = 0L))

#'@importFrom httr GET accept_json content
#'@importFrom stringr str_c
#'
hpo_api_get <- function(extension,
                        base_url = "https://hpo.jax.org/api/hpo/")
{
  url <- str_c(base_url, extension)
  
  response <- GET(url, accept_json())
  
  if (response$status_code == 200) {
    return(content(response))
  }
  warning('GET ', url, ' failed.')
  return(NULL)
}

#'@importFrom dplyr tibble as_tibble bind_rows "%>%"
#'@importFrom urltools url_encode
#'@importFrom memoise memoise 
#'@importFrom purrr map_df
#'
get_hpo_term <- function(hpo_id)
{
  assert_that(is.character(hpo_id),
              all(str_starts(hpo_id, 'HP:'), na.rm = TRUE))
  
  mapper <- cavalier_cache$hpo_api_get_term_mapper
  
  if (is.null(mapper)) {
    mapper <- memoise(
      function(extension) {
        if (!is.na(extension)) {
          result <- hpo_api_get(extension)
          if (!is.null(result)) {
            data <-
              result$details %>%
              map(~ `if`(is_list(.), list(unlist(.)), .)) %>% 
              as_tibble() %>% 
              mutate(parents = list(bind_rows(result$relations$parents)),
                     children = list(bind_rows(result$relations$children)))
            return(data)
          } 
        }
        return(tibble(name = NA_character_))
      })
    cavalier_cache$hpo_api_get_term_mapper <- mapper
  }
  
  map_df(str_c('term/', url_encode(hpo_id)), mapper)
}

#'@importFrom rlang is_integerish 
get_hpo_gene <- function(entrez_id)
{
  assert_that(is_integerish(entrez_id))
  
  mapper <- cavalier_cache$get_hpo_gene_mapper
  
  if (is.null(mapper)) {
    mapper <- memoise(
      function(entrez_id) {
        if (!is.na(entrez_id)) {
          result <- hpo_api_get(str_c('gene/', entrez_id))
          if (!is.null(result)) {
            data <- 
              as_tibble(result$gene) %>% 
              mutate(hpo_terms = list(bind_rows(result$termAssoc)),
                     diseases = list(bind_rows(result$diseaseAssoc)))
            return(data)
          } 
        }
        return(tibble(entrezGeneId = entrez_id))
      })
    cavalier_cache$get_hpo_gene_mapper <- mapper
  }
  
  map_df(entrez_id, mapper)
}

get_hpo_disease <- function(disease_id)
{
  assert_that(
    is.character(disease_id),
    all(str_starts(disease_id, '(OMIM)|(ORPHA):'), na.rm = TRUE))
  
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
                hpo_terms = list(map_df(result$catTermsMap, function(data) {
                  tibble(catLabel = data$catLabel,
                         terms = list(bind_rows(data$terms))) %>% 
                    unnest(terms)
                })),
                genes = list(bind_rows(result$geneAssoc)))
            return(data)
          } 
        }
        return(tibble(diseaseId = disease_id))
      })
    cavalier_cache$get_hpo_disease_mapper <- mapper
  }
  
  map_df(disease_id, mapper)
}