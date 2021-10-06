
#' @importFrom dplyr row_number mutate
#' @importFrom purrr map2
#' @importFrom readr write_tsv
#' @export
get_web_list_version <- function(id, save = NULL)
{
  assert_that(is_character(id))
  
  result <-
    tibble(id = id) %>% 
    mutate(row_num = row_number(),
           source = str_extract(id, '^[^:]+')) %>% 
    nest(data = -source) %>% 
    mutate(data = map2(source, data, function(source, data) {
      if (source %in% c('PAA', 'PAE')) {
        get_panelapp_panels(source) %>% 
          select(id, version) %>% 
          left_join(data, ., 'id')
      } else if (source == 'HP') {
        data %>% 
          mutate(version = latest_hpo_build_num())
      } else {
        data %>% 
          mutate(version = NA_character_)
      }
    })) %>% 
    unnest(data) %>% 
    arrange(row_num) %>% 
    select(id, version)
  
  if (!is.null(save)) {
    write_tsv(result, save)
    invisible(result)
  } else {
    result
  }
}


#' @export
get_web_list <- function(id, save = NULL, secure = TRUE)
{
  assert_that(is_scalar_character(id))
  
  if (!secure) { insecure() }
  
  result <- 
    switch(str_extract(id, '^[^:]+'),
           HP = get_hpo_gene_list(id),
           PAA = get_panelapp_gene_list(id),
           PAE = get_panelapp_gene_list(id))
  
  if (!is.null(save)) {
    write_tsv(result, save)
    invisible(result)
  } else {
    result
  }
}
