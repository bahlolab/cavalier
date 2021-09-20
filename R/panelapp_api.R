
get_panelapp_url <- function(source = c('PAA', 'PAE')) 
{
  source <- match.arg(source)
  if (source == 'PAA') {
    return("https://panelapp.agha.umccr.org/")
  }
  if (source == 'PAE') {
    return("https://panelapp.genomicsengland.co.uk/")
  }
}

#'@importFrom dplyr tibble as_tibble bind_rows mutate bind_cols
#'@importFrom tidyr unnest
#'@importFrom httr GET accept_json content
#'@importFrom purrr map_df map
#'@importFrom stringr str_c
#'@importFrom rlang is_scalar_integerish
#'@export
get_panelapp_panel <- function(id, min_confidence = 2L) 
{
  assert_that(is_scalar_character(id),
              !is.na(id),
              str_detect(id, '^PA[AE]:\\d+'))
  
  url <- str_c( get_panelapp_url(str_extract(id, '^PA[AE]')), 
                'api/v1/panels/', str_extract(id, '\\d+'), '/')
  
  GET(url, accept_json()) %>%
    content() %>%
    (function(x) {
      tibble(id = x$id,
             names = x$name,
             version = x$version,
             gene_data = list(map_df(x$genes, function(y) {
               cols1 <- c("entity_type", "entity_name", "confidence_level",
                          "mode_of_pathogenicity", "mode_of_inheritance")
               cols2 <- c('biotype', 'hgnc_id', 'gene_name', 'gene_symbol',
                          'hgnc_symbol', 'hgnc_release', 'hgnc_date_symbol_changed')
               as_tibble(map(y[cols1], ~ `if`(is.null(.), NA_character_, .))) %>%
                 bind_cols(
                   as_tibble(map(y$gene_data[cols2], ~ `if`(is.null(.), NA_character_, .)))) %>%
                 mutate(alias = list(y$gene_data$alias),
                        phenotypes = list(unlist(y$phenotypes)),
                        evidence = list(unlist(y$evidence))) %>% 
                 mutate(panel_id = y$panel$id,
                        panel_name =y$panel$name,
                        panel_version = y$panel$version)
             }))) %>% 
        (function(x) {
          if (!'panel_id' %in% names(x)) {
            x$panel_id <- x$id
          }
          if (!'panel_name' %in% names(x)) {
            x$panel_name <- x$names
          }
          if (!'panel_version' %in% names(x)) {
            x$panel_version <- x$version
          }
          x
        })
    }) %>%
    unnest(gene_data) %>% 
    filter(entity_type == 'gene',
           confidence_level >= min_confidence) %>% 
    mutate(status = case_when(confidence_level == 3 ~ 'GREEN',
                              confidence_level == 2 ~ 'AMBER',
                              confidence_level == 1 ~ 'RED'),
           inheritance = str_extract(mode_of_inheritance, 
                                     '^(BIALLELIC)|(MONOALLELIC)|(BOTH)|(X-LINKED)')) %>% 
    select(list_id = panel_id,
           list_name = panel_name,
           version = panel_version,
           gene = entity_name,
           inheritance,
           status)
}