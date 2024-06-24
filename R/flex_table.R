
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
flex_table_trans <- function(data, url_df = NULL, ncol = 1L)
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
    (function(x) {
      split <- parallel::splitIndices(nrow(x), ncol)
      nr <- max(lengths(split))
      map2(split, seq_along(split), function(ix, i) {
        pad_df(x[ix, ], nr) %>% rename_with(~ str_c(.,'.', i)) 
      }) %>% do.call(bind_cols, .)
    }) %>% 
    flextable(col_keys = flatten_chr(
      map(seq_len(ncol), ~ str_c(c('name', 'value'), '.', .)))) %>% 
    (function(x) {
      `if`(has_url,
           seq_len(ncol) %>% 
             reduce(function(ft, i) {
               url_col <- str_c('url.', i)
               text_col <- str_c('value.', i)
               compose(ft, j = text_col,
                       value = as_paragraph(
                         hyperlink_text(x = !!sym(text_col),
                                        url = !!sym(url_col)
                         )))
             },
             .init = x),
           x)
    }) %>% 
    delete_part(part = "header") %>% 
    theme_zebra(even_header = 'white', even_body = 'white') %>% 
    italic(j = (seq_len(ncol) * 2) - 1) %>% 
    bold(j = (seq_len(ncol) * 2) - 1) %>% 
    colformat_char(j = (seq_len(ncol) * 2) - 1, suffix = ':') %>% 
    align(j = (seq_len(ncol) * 2) - 1, align = 'right', part = 'all')  %>%
    flextable::vline(j = seq_len(ncol-1) * 2, border = fp_border_default(color = '#CFCFCF'))
}


#' @importFrom flextable fontsize autofit dim_pretty width height_all fit_to_width
fit_flex_table <- function(ft, height, width,
                           font_size = 12L,
                           min_font_size = 8L,
                           add_h = -0.05,
                           add_w = -0.1,
                           max_lines = 10L,
                           expand_cols = TRUE,
                           expand_rows = FALSE)
{
  while (font_size >= min_font_size) {
    
    ft <- fontsize(ft, size = font_size, part = 'all') 
    cell_dim <- 
      cell_dims_wrapped(ft, max_lines = max_lines) %>% 
      mutate(part = ordered(part, c('header', 'body', 'footer')),
             col_id = ordered(col_id, ft$col_keys))
    
    width_constraint <-
      cell_dim %>% 
      mutate(width_min = map_dbl(dims, ~ min(.$width))) %>% 
      group_by(col_id) %>% 
      summarise(width_min = max(width_min) + add_w, .groups = 'drop') %>% 
      mutate(width_max = width - sum(width_min) + width_min)
    
    if (sum(width_constraint$width_min) > width & font_size > min_font_size) {
      # too wide, try smaller font
      font_size <- font_size - 1L
      next()
    }
    
    # drop cell dims based on min_width and width_max
    cell_dim <-
      cell_dim %>% 
      left_join(width_constraint, by = 'col_id') %>% 
      mutate(dims_f = pmap(., function(width_min, width_max, dims, ...) {
        dims %>% 
          mutate(width = width + add_w) %>% 
          filter(width >= width_min | width == suppressWarnings(max(width[width < width_min])),
                 width <= width_max)
      })) %>% 
      select(part, row_id, col_id, dims, dims_f)
    
    if (any(map_int(cell_dim$dims_f, nrow) == 0)) {
      if (font_size > min_font_size) {
        font_size <- font_size - 1L
        next()
      } 
      cell_dim <-
        mutate(cell_dim, dims = if_else(
          map_int(dims_f, nrow) == 0,
          map(dims, ~ slice(arrange(., width * height), 1)),
          dims_f
        ))
    }
    
    height_constraint <-
      cell_dim %>% 
      mutate(height_min = map_dbl(dims, ~ min(.$height))) %>% 
      group_by(part, row_id) %>% 
      summarise(height_min = max(height_min) + add_h, .groups = 'drop') %>% 
      mutate(height_max = height - sum(height_min) + height_min)
    
    if (sum(height_constraint$height_min) > height & font_size > min_font_size) {
      # too tall, try smaller font
      font_size <- font_size - 1L
      next()
    }
    
    # drop cell dims based on  min and max height
    cell_dim <-
      cell_dim %>% 
      left_join(height_constraint, by = c('part',  'row_id')) %>% 
      mutate(dims_f = pmap(., function(height_min, height_max, dims, ...) {
        dims %>% 
          mutate(height = height + add_h) %>% 
          filter(height >= height_min | height == suppressWarnings(max(height[height < height_min])),
                 height <= height_max)
      })) %>% 
      select(part, row_id, col_id, dims, dims_f)
    
    if (any(map_int(cell_dim$dims_f, nrow) == 0)) {
      if (font_size > min_font_size) {
        font_size <- font_size - 1L
        next()
      }
      cell_dim <-
        mutate(cell_dim, dims = if_else(
          map_int(dims_f, nrow) == 0,
          map(dims, ~ slice(arrange(., width * height), 1)),
          dims_f
        ))
    }
    
    # try all combinations of cells with multiple dims
    dimopts <-
      cell_dim %>% 
      split.data.frame(seq_along(.$part)) %>% 
      map(~ unnest(., dims)) %>% 
      do.call(expand_grid, .) %>% 
      as.list() %>% 
      map(mutate, id = seq_along(part)) %>% 
      bind_rows() %>% 
      arrange(id, part, row_id, col_id) %>% 
      (function(x) full_join(
        group_by(x, id, part, row_id) %>% 
          summarise(height = max(height),
                    .groups = 'drop') %>% 
          select(id, part, heights = height) %>% 
          chop(c(part, heights)),
        group_by(x, id, col_id) %>% 
          summarise(width = max(width),
                    .groups = 'drop') %>% 
          select(id, widths = width) %>% 
          chop(widths),
        by = 'id'
      )) %>% 
      mutate(tot_height = map_dbl(heights, sum),
             tot_width = map_dbl(widths, sum),
             prop_height = tot_height / height,
             prop_width =  tot_width / width)
    
    
    if (!any(dimopts$prop_height <= 1 & dimopts$prop_width <= 1) & 
        font_size > min_font_size) {
      # no viable dims
      font_size <- font_size - 1L
      next()
    }
    
    if (any(dimopts$prop_height <= 1 & dimopts$prop_width <= 1)) {
      dimopts <-  
        dimopts %>% 
        filter(dimopts$prop_height <= 1,
               dimopts$prop_width <= 1)
    }
    # choose option using least area
    opt_dim <-
      dimopts %>% 
      arrange(prop_height * prop_width) %>% 
      slice(1) %>% 
      select(part, heights, widths, tot_height, tot_width)
    
    break()
  }
  
  
  if (expand_rows || opt_dim$tot_height > height) {
    opt_dim$heights[[1]] <-
      height * opt_dim$heights[[1]] / sum(opt_dim$heights[[1]])
  }
  if (expand_cols || opt_dim$tot_width > width) {
    opt_dim$widths[[1]] <-
      width * opt_dim$widths[[1]] / sum(opt_dim$widths[[1]])
  }
  
  opt_dim %>% 
    select(part, heights) %>% 
    unnest(c(part, heights)) %>% 
    mutate(part = as.character(part)) %>% 
    chop(heights) %>% 
    pwalk(function(part, heights) {
      ft[[part]]$rowheights <<- heights
    })
  
  c('header', 'body', 'footer') %>% 
    walk(function(part) {
      ft[[part]]$colwidths <<- opt_dim$widths[[1]]
    })
  
  return(ft)
}

# return heights and widths for flextable with wrapping
#' @importFrom tidyr nest
cell_dims_wrapped <- function(ft,
                              max_lines = Inf)
{
  map_df(c("header", "body", "footer"), function(part_name) {
    
    part <- ft[[part_name]]
    
    if (length(unlist(part$spans)) == 0) {
      return(tibble())
    }
    
    non_txt_dim <- 
      list(flextable:::dim_cells(part),
           flextable:::dim_paragraphs(part)) %>% 
      map(function(x) map(x, as_tibble) %>% 
            bind_cols() %>% 
            mutate(row_id = row_number())) %>% 
      bind_rows() %>% 
      pivot_longer(-row_id,
                   names_to = c('.value', 'col_id'),
                   names_pattern = '([^.]+)\\.(.+)') %>% 
      group_by(row_id, col_id) %>% 
      summarise(across(everything(), sum),
                .groups = 'drop')
    
    flextable:::fortify_content(part$content, default_chunk_fmt = part$styles$text) %>% 
      (function(x) {
        tibble(row_id = x$.row_id,
               col_id = as.character(x$.col_id),
               dims =  select(x, txt, font.size, font.family, bold, italic) %>% 
                 pmap(., txt_dim_wrap, max_lines = max_lines))
      }) %>% 
      unnest(dims) %>% 
      left_join(non_txt_dim, by = c("row_id", "col_id"), suffix = c('', '_nt')) %>% 
      mutate(width = width + width_nt,
             height = height + height_nt) %>% 
      mutate(part = part_name) %>% 
      select(part, row_id, col_id, height, width) %>% 
      nest(dims = c(height, width))
  })
  
}

txt_dim_wrap <- function(txt, font.size, font.family, bold, italic,
                         max_lines = Inf) {
  
  if (is.na(txt) || nchar(txt) == 0) {
    return(tibble(width=0, height=0))
  }
  
  tibble(nchar = seq.int(nchar(txt))) %>% 
    mutate(text = map2_chr(txt, nchar, ~ str_wrap(.x, width = .y))) %>% 
    select(text) %>% 
    distinct() %>% 
    filter(str_count(text, '\\n') <= max_lines) %>% 
    mutate(id = seq_along(text)) %>% 
    separate_rows(text, sep = '\n') %>% 
    chop(id) %>% 
    (function(x) {
      gdtools::m_str_extents(x$text, 
                             fontname = font.family,
                             fontsize = font.size, 
                             bold = bold,
                             italic = italic) %>% 
        set_colnames(c('width', 'height')) %>%
        as_tibble() %>% 
        bind_cols(x, .) %>% 
        mutate(across(c(width, height), ~ . / 72))
    }) %>% 
    select(-text) %>% 
    unchop(id) %>% 
    group_by(id) %>% 
    summarise(n_lines = n(),
              height = sum(height),
              width = max(width),
              .groups = 'drop') %>% 
    arrange(n_lines, width) %>% 
    group_by(n_lines) %>% 
    slice(1) %>% 
    ungroup() %>% 
    select(width, height) %>% 
    arrange(width, desc(height)) %>% 
    # ensure monotonic
    filter(map_lgl(seq_along(height), function(i) 
      all(width[i] < width[-i] | ( width[i] >= width[-i] & height[i] <= height[-i]))))
}
