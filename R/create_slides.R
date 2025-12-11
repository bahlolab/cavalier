#' @importFrom officer read_pptx add_slide ph_with ph_location_type ph_location_template external_img 
#' @importFrom tibble rownames_to_column
#' @export
create_slides <- function(
    slide_layout,
    slide_data,
    output = 'cavalier_slides.pptx',
    slide_template = get_slide_template()
)
{
  # check args
  assert_that(
    is.data.frame(slide_layout),
    is.data.frame(slide_data),
    is_scalar_character(output),
    is_scalar_character(slide_template) && file.exists(slide_template)
  )
  
  SLIDES <- read_pptx(slide_template)
  n_slides <- n_distinct(slide_layout$slide_num)
  
  for (i in seq_len(nrow(slide_data))) {
    
    DATA <-  map(slide_data[i,], function(x) { if(is.list(x) && length(x) == 1) { x[[1]]} else { x} })
    message('Adding slides ', i, ' of ', nrow(slide_data), ': "', DATA$TITLE, '"')
    
    for (j in seq_len(n_slides)) {
      
      TITLE <- DATA$TITLE
      if (n_slides > 1) {
        TITLE <- str_c(TITLE, ' (', j, '/', n_slides, ')')
      }
      
      SLIDES <-
        SLIDES %>% 
        add_slide(layout = "Title and Content") %>% 
        ph_with(value = TITLE, location = ph_location_type(type = "title"))
      
      slide_layout %>%
        filter(slide_num == j) %>%
        pwalk(function(element, x_left, width, y_top, height, transpose, ...) {
          value <- DATA[[element]]
          
          if (is.data.frame(value)) {
            value <-
              value %>% 
              flex_table(transpose = transpose) %>% 
              fit_flex_table(width = width, height = height, expand_rows = transpose)
          } else if (is(value, 'flextable')) {
            value <- fit_flex_table(value, width = width, height = height, expand_rows = FALSE)
          }
          # add item to slides
          if (!is.null(value)) {
            SLIDES <-
              SLIDES %>% 
              ph_with(value = value,
                      location = ph_location_template(
                        left = x_left,
                        top = y_top,
                        width = width,
                        height = height))
          }
        })
    }
  }
  
  print(SLIDES, target = output)
  
  # no longer required?
  # re_encode_pptx_hlinks(output)
  
  message("Created slides: ", output)
  
  return(invisible(output))
}

#' @export
get_slide_template <- function() {
  system.file("ppt", "template.pptx", package = "cavalier")
}

#' @importFrom rlang dots_list is_scalar_double
#' @export
slide_layout <- function(...,
                         heights = NULL,
                         title_height = 0.1,
                         slide_height = 7.5,
                         slide_width = 13.333,
                         pad = 0.02,
                         transpose = character(),
                         slide_num = 1L)
{
  rows <- dots_list(...)
  
  assert_that(
    is_scalar_double(title_height), title_height >= 0, title_height < 1,
    is_scalar_double(slide_height),
    is_scalar_double(slide_width),
    is_scalar_double(pad), pad >= 0, pad < 1,
    length(rows) > 0,
    # all(map_lgl(rows, is_valid_row)),
    is.null(heights) | (is_number(heights) & length(heights) == length(rows)))
  
  
  # get coordinates for each element
  layout_df <-
    map_df(rows, function(row) {
      `if`(is_named(row), 
           tibble(element = names(row),
                  width = row / sum(row)),
           tibble(element = row,
                  width = rep(1/length(row), length(row)))) %>% 
        pad_x(pad = pad * slide_width, 
              slide_width = slide_width) %>% 
        nest(data = everything())
    }) %>% 
    mutate(height = `if`(is.null(heights),
                         rep(1/length(rows), length(rows)),
                         heights / sum(heights))) %>% 
    pad_y(pad = pad * slide_width, 
          title_height = title_height * slide_height,
          slide_height = slide_height) %>% 
    unnest(data) %>% 
    mutate(slide_num = slide_num,
           transpose = element %in% transpose
    )
  
  return(layout_df)
}

pad_x <- function(col_df, pad, slide_width) {
  mutate(col_df, 
         abs = FALSE,
         row_num = row_number() * 2L) %>% 
    bind_rows(tibble(element = 'pad',
                     width = pad,
                     abs = TRUE,
                     row_num = seq.int(1, nrow(col_df) * 2 +1, by = 2))) %>% 
    arrange(row_num) %>% 
    mutate(width = if_else(abs, width, width * (slide_width - sum(width[abs]))),
           x_right = cumsum(width),
           x_left = x_right - width,
           x_cen = (x_right + x_left) / 2) %>% 
    select(element, width, x_left, x_cen, x_right) %>% 
    filter(element != 'pad')
}

pad_y <- function(row_df, pad, title_height, slide_height) {
  mutate(row_df, 
         abs = FALSE,
         row_num = row_number() * 2L) %>% 
    bind_rows(tibble(data = list(NULL),
                     height = pad,
                     abs = TRUE,
                     row_num = seq.int(1, nrow(row_df) * 2 +1, by = 2))) %>% 
    add_row(row_num = 0L,
            data = NULL, 
            abs = TRUE,
            height = title_height) %>% 
    arrange(row_num) %>% 
    mutate(height = if_else(abs, height, height * (slide_height - sum(height[abs]))),
           y_bot = cumsum(height),
           y_top = y_bot - height,
           y_cen = (y_bot + y_top) / 2) %>%
    select(data, height, y_bot, y_cen, y_top) %>%
    filter(map_lgl(data, ~ ! is.null(.)))
}

# # Resolved in recent Officer versions?
# re_encode_pptx_hlinks <- function(target) 
# {
#   hlink_type <- "http://schemas.openxmlformats.org/officeDocument/2006/relationships/hyperlink"
#   tmp_dir <- tempfile(pattern = '.', tmpdir = '.')
#   zip::unzip(target, exdir = tmp_dir)
#   
#   changed <- FALSE
#   
#   list.files(file.path(tmp_dir, 'ppt/slides/_rels'),
#              full.names = TRUE) %>% 
#     walk(function(file) {
#       content <- 
#         xml2::read_xml(file) %>% 
#         xml2::as_list()
#       
#       rep_rel <-
#         map(content$Relationships, function(item) {
#           atrib <- attributes(item)
#           if (atrib$Type == hlink_type) {
#             if (str_detect(atrib$Target, '%')) {
#               atrib$Target <- 
#                 utils::URLdecode(atrib$Target) %>% 
#                 utils::URLencode()
#               changed <<- TRUE
#             }
#           }
#           attributes(item) <- atrib
#           item
#         })
#       
#       attributes(rep_rel) <- attributes(content$Relationships)
#       content$Relationships <- rep_rel
#       
#       if (changed) {
#         xml2::write_xml(xml2::as_xml_document(content), file)
#       }
#     })
#   
#   if (changed) {
#     officer::pack_folder(tmp_dir, target)
#   }
#   unlink(tmp_dir, recursive = TRUE)
#   
#   invisible(NULL)
# }

