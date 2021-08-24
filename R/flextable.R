#' @importFrom flextable flextable delete_part theme_zebra italic bold colformat_char 
#' @importFrom flextable align autofit fit_to_width add_header_lines
flextable_reg <- function(table_data,
                          title = NULL,
                          digits = 3) {
  
  table_data %>% 
    mutate(across(where(is.double), ~ as.character(round(., digits))),
           across(everything(), as.character),
           across(everything(), ~replace_na(., 'NA'))) %>% 
    flextable() %>% 
    add_header_lines(title) %>% 
    theme_zebra(even_header = 'white', even_body = 'white') %>% 
    align(i =1, part = 'header', align = 'center')
  
}

#' @importFrom flextable flextable delete_part theme_zebra italic bold colformat_char align autofit fit_to_width
flextable_trans <- function(table_data,
                            digits = 3) {
  
  table_data %>% 
    mutate(across(where(is.double), ~ as.character(round(., digits))),
           across(everything(), as.character),
           across(everything(), ~replace_na(., 'NA'))) %>% 
    pivot_longer(everything()) %>% 
    flextable() %>% 
    delete_part(part = "header") %>% 
    theme_zebra(even_header = 'white', even_body = 'white') %>% 
    italic(j = 1) %>% 
    bold(j = 1) %>% 
    colformat_char(j = 1, suffix = ':') %>% 
    align(j = 1, align = 'right', part = 'all')
}

#' @importFrom flextable fontsize autofit dim_pretty width height_all fit_to_width
flextable_fit <- function(ft, width, height,
                          start_size = 11,
                          min_size = 5,
                          max_size = 16,
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
