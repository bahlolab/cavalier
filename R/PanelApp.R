
#' Get URL for PanelApp API
get_panelapp_url <- function(source) 
{
  panelapp_urls <- get_cavalier_opt("panelapp_urls")
  
  assert_that(
    is_scalar_character(source),
    is.list(panelapp_urls),
    source %in% names(panelapp_urls),
    is_scalar_character(panelapp_urls[[source]])
  )
  
  return(panelapp_urls[[source]])
}

#' Get table of all current PanelApp panels
#'@importFrom dplyr tibble as_tibble bind_rows mutate bind_cols
#'@importFrom tidyr unnest
#'@importFrom httr GET accept_json content
#'@importFrom purrr map_df map
#'@importFrom stringr str_c
#'@importFrom rlang is_scalar_integerish
#'@export
get_panelapp_panels <- function(source) 
{
  url <- str_c(get_panelapp_url(source), 'api/v1/panels')
  
  res <- tibble()
  while (TRUE) {
    response <- retry('GET', url, accept_json()) %>% content()
    results <-
      response$results %>%
      map_df(function(x) {
        if (is.null(x$hash_id)) {
          x$hash_id <- NA_character_
        }
        types <- x$types
        x$types <- NULL
        stats <- x$stats
        x$stats <- NULL
        if (length(x$relevant_disorders) == 0) {
          x$relevant_disorders <- list(character())
        } else {
          x$relevant_disorders <- list(x$relevant_disorders)
        }
        x <- c(x,
               setNames(stats, str_c('stats_', names(stats))),
               list(types = list(map_df(types, as_tibble))))
        as_tibble(x)
      })
    
    res <- bind_rows(res, results)
    if (is.null(response$`next`)) { break }
    url <- response$`next`
  }
  
  mutate(res, id = str_c(source, ':', id))
  
}

#' Get latest version of PanelApp Panel (from server or cache)
get_panelapp_panel_version <- function(
    id,
    db_mode = get_cavalier_opt("database_mode")
) 
{
  
  
  assert_that(is_scalar_character(id),
              !is.na(id),
              str_detect(id, '^[A-z]+:\\d+'))
  
  func_online <- function() {
    source <- str_extract(id, '^[A-z]+')
    get_panelapp_panels(source = source) %>% 
      filter(id == !!id) %>% 
      pull(version) %>% 
      first()
  }
  
  get_version(
    resource_name  = id,
    cache_name = str_replace(id, ':', '_'),
    cache_subdir = 'PanelApp',
    func_online = func_online,
    db_mode = db_mode
  )
}

#' Get gene panel from PanelApp API
#'@importFrom dplyr tibble as_tibble bind_rows mutate bind_cols
#'@importFrom tidyr unnest
#'@importFrom httr GET accept_json content RETRY
#'@importFrom purrr map_df map
#'@importFrom stringr str_c
#'@importFrom rlang is_scalar_integerish
#'@export
get_panelapp_panel <- function(id, version = NULL) 
{
  assert_that(is_scalar_character(id),
              !is.na(id),
              str_detect(id, '^[A-z]+:\\d+'),
              is.null(version) | is_scalar_character(version))
  
  source <- str_extract(id, '^[A-z]+')
  
  if (is.null(version)) {
    version <- get_panelapp_panel_version(id)
  }
  
  url <- str_c(
    get_panelapp_url(source), 
    'api/v1/panels/', str_extract(id, '\\d+'), '/',
    `if`(!is.null(version), str_c('?version=', version))
  )
  
  fun <- function() {
    retry('GET', url, accept_json()) %>%
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
                          panel_name = y$panel$name,
                          panel_version = y$panel$version)
               }))) 
      }) %>%
      unnest(gene_data) %>% 
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
      }) %>% 
      bind_rows(tibble(
        entity_type = character(), 
        entity_name = character(), 
        hgnc_id = character(),
        hgnc_symbol = character(),
        gene_symbol = character(),
        confidence_level = character(),
        mode_of_inheritance = character(),
      )) %>% 
      mutate(confidence_level = as.integer(confidence_level)) %>% 
      filter(entity_type == 'gene') %>% 
      mutate(panel_id = str_c(source, ':', panel_id),
             status = case_when(confidence_level == 3 ~ 'GREEN',
                                confidence_level == 2 ~ 'AMBER',
                                confidence_level == 1 ~ 'RED'),
             inheritance = str_extract(mode_of_inheritance, 
                                       '^(BIALLELIC)|(MONOALLELIC)|(BOTH)|(X-LINKED)')) %>% 
      mutate(inheritance = case_when(
        inheritance == 'BIALLELIC' ~ 'AR',
        inheritance == 'MONOALLELIC' ~ 'AD',
        inheritance == 'BOTH' ~ 'AR/AD',
        inheritance == 'X-LINKED' ~ 'XL')) %>% 
      select(list_id = panel_id,
             list_name = panel_name,
             version = panel_version,
             hgnc_id,
             gene = hgnc_symbol,
             panelapp_symbol = gene_symbol,
             inheritance,
             status,
             confidence_level)
  }
  
  cache(
    fun = fun,
    name = str_replace(id, ':', '_'),
    version = version,
    subdir = 'PanelApp'
  )
}

#' Get version history for PanelApp panel
get_panelapp_version_history <- function(id) 
{
  assert_that(is_scalar_character(id),
              !is.na(id),
              str_detect(id, '^[A-z]+:\\d+'))
  
  source <- str_extract(id, '^[A-z]+')
  
  url <- str_c(
    get_panelapp_url(source), 
    'api/v1/panels/',
    str_extract(id, '\\d+'),
    '/activities'
  )
  
  response <- 
    retry('GET', url, accept_json()) %>%
    content()
  
  versions <- unique(map_chr(response, 'panel_version'))
  
  return(versions)
}

#' Store local cache of all PanelApp panels
build_panelapp_cache <- function(sources = names(get_cavalier_opt('panelapp_urls'))) {
  
  for(source in sources) {
    get_panelapp_panels(source = source) %>% 
      select(id, version) %>% 
      # head(10) %>% # testing only
      pwalk(function(id, version)
        invisible(get_panelapp_panel(id=id, version=version))
      )
  }
}
