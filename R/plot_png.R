#' @importFrom png readPNG writePNG
#' @importFrom stringr str_c str_remove
#' @importFrom abind abind
crop_png <- function(
    input_png,
    output_png,
    crop_left = 0, 
    crop_right = 0,
    crop_top = 0,
    crop_bot = 0
) {
  assert_that(
    is_character(input_png),
    all(file.exists(input_png)),
    is_character(output_png),
    length(input_png) == length(output_png),
    is_scalar_integerish(crop_left),
    is_scalar_integerish(crop_right),
    is_scalar_integerish(crop_top),
    is_scalar_integerish(crop_bot)
  )
  
  map2_df(input_png, output_png, function(input, output) {
    png <- readPNG(input)
    # crop left, right and top and bottom
    png <- png[seq.int(crop_top  +1, dim(png)[1] - crop_bot),
               seq.int(crop_left +1, dim(png)[2] - crop_right),]
    writePNG(png, target = output)
    tibble( width = dim(png)[2], height = dim(png)[1], png = output)
  })
}


#' @importFrom ggplot2 ggplot theme_void coord_fixed aes scale_x_continuous scale_y_continuous
#' @importFrom ggimg geom_rect_img
#' @export
plot_png_facets <- function(
    id_png_tbl,
    crop_left = 0, 
    crop_right = 0,
    max_cols = 3,
    base_size = 14
) {
  
  if (crop_left + crop_right > 0) {
    id_png_tbl <-
      id_png_tbl %>% 
      mutate(cropped = str_c(png, '.crop.png'))
    dims <-
      id_png_tbl %>% 
      with(crop_png(png, cropped, crop_left, crop_right)) %>% 
      select(width, height) %>% 
      map(mean) %>% 
      map(round)
    
    id_png_tbl <-
      id_png_tbl %>% 
      select(id, png = cropped)
  } else{
    dims <-
      map_df(id_png_tbl$png, function(x) {
        png <- readPNG(x)
        tibble( width = dim(png)[2], height = dim(png)[1])
      }) %>% 
      map(mean) %>% 
      map(round)
  }

  p <-
    id_png_tbl %>% 
    ggplot() + 
    geom_rect_img(aes(xmin = 0, xmax = dims$width, ymin = 0, ymax = dims$height, img = png)) +
    facet_wrap(~id, strip.position = 'top',
               ncol = min(nrow(id), max_cols)) + 
    theme_void(base_size = base_size) + 
    coord_fixed() +
    scale_x_continuous(limits = c(0, dims$width), expand = c(0.01, 0.01)) +
    scale_y_continuous(limits = c(0, dims$height), expand = c(0.01, 0.01)) 
  return(p)
}
