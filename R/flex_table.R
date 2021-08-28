
#' @importFrom flextable flextable delete_part theme_zebra italic bold colformat_char fit_to_width
#' @importFrom flextable align autofit add_header_lines border_inner fp_border_default
#' @importFrom flextable compose as_paragraph hyperlink_text
flex_table <- function(data,
                       transpose = FALSE,
                       round = TRUE,
                       digits = 3)
{
  assert_that(is.data.frame(data),
              is_bool(transpose),
              is_bool(round),
              is_scalar_integerish(digits))
  url_df <-
    tibble(url_col = keep(names(data), ~ str_ends(., '_url')),
           text_col = str_remove(url_col, '_url')) %>% 
    filter(text_col %in% names(data))
  
  data <- 
    data %>% 
    (function(x) 
      `if`(round,
           mutate(x, across(where(is.double), ~ as.character(round(., digits)))),
           x) 
    ) %>% 
    mutate(across(everything(), as.character),
           across(!all_of(url_df$url_col), ~replace_na(., 'NA')))
  
  if (transpose) {
    return(flex_table_trans(data, url_df))
  }
  
  data %>% 
    flextable(col_keys = setdiff(names(data), url_df$url_col)) %>% 
    (function(ft) {
      reduce(seq_len(nrow(url_df)),
             function(ft, i) {
               url_col <- url_df$url_col[i]
               text_col <- url_df$text_col[i]
               compose(ft, j = text_col,
                       value = as_paragraph(
                         hyperlink_text(x = !!sym(text_col),
                                        url = !!sym(url_col)
                         )))
             },
             .init = ft)
    }) %>% 
    theme_zebra(even_body = 'white') %>% 
    border_inner(border = fp_border_default(color = '#CFCFCF'))
}

#' @importFrom flextable flextable delete_part theme_zebra italic bold colformat_char 
#' @importFrom flextable align autofit fit_to_width
flex_table_trans <- function(data, url_df = NULL)
{
  assert_that(is.data.frame(data),
              nrow(data) <= 1,
              is.null(url_df) | is.data.frame(url_df))
  
  has_url <- !is.null(url_df) & nrow(url_df)
  
  data %>% 
    select(!all_of(url_df$url_col)) %>% 
    pivot_longer(everything()) %>% 
    (function(x) {
      `if`(has_url,
           data %>% 
             select(all_of(url_df$url_col)) %>% 
             pivot_longer(everything(),  values_to = 'url') %>% 
             mutate(name = str_remove(name, '_url$')) %>% 
             left_join(x, ., by = 'name'),
           x)
    }) %>% 
    flextable(col_keys = c('name', 'value')) %>% 
    (function(x) {
      `if`(has_url,
           compose(x, j = 2,
                   value = as_paragraph(
                     hyperlink_text(x = value,
                                    url = url
                     ))),
           x)
    }) %>% 
    delete_part(part = "header") %>% 
    theme_zebra(even_header = 'white', even_body = 'white') %>% 
    italic(j = 1) %>% 
    bold(j = 1) %>% 
    colformat_char(j = 1, suffix = ':') %>% 
    align(j = 1, align = 'right', part = 'all')
}

#' @importFrom flextable fontsize autofit dim_pretty width height_all fit_to_width
fit_flex_table <- function(ft, width, height,
                          start_size = 12,
                          min_size = 5,
                          max_size = 14,
                          max_row_height = 0.33) {
  # goal is to shrink until both width and height are less than dim_pretty
  curr_size <- start_size
  ft <- 
    fontsize(ft, size = curr_size, part = 'all') %>% 
    autofit()    
  dims <- dim_pretty(ft) %>% map(sum)
  # too small 
  while (curr_size < max_size & dims$heights < height & dims$widths < width) {
    curr_size <- curr_size + 1L
    ft <- 
      fontsize(ft, size = curr_size, part = 'all') %>%
      autofit()
    dims <- dim_pretty(ft) %>% map(sum)
  }
  # too big 
  while (curr_size > min_size & (dims$heights > height | dims$widths > width)) {
    curr_size <- curr_size - 1L
    ft <- 
      fontsize(ft, size = curr_size, part = 'all') %>% 
      autofit()
    dims <- dim_pretty(ft) %>% map(sum)
  }
  dims <- dim_pretty(ft)
  dims$widths <- (dims$widths / sum(dims$widths)) * width
  height <- {height / length(dims$heights) } %>% min( max_row_height)
  ft %>%
    width(seq_along(dims$widths), dims$widths) %>% 
    height_all(height)
}
