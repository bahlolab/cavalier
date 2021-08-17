# Return directory name that ends in "/" (for pasting strings)
endslash_dirname <- function(dirname)
{
    if (!endsWith(dirname, "/")) {
        dirname <- paste0(dirname, "/")
    }
    return(dirname)
}

# Insert a newline every N characters (e.g., for long genome changes)
newline_every_n_chars <- function(x, n)
{
    if (nchar(x) > n) {
        x0 <- substr(x, 1, n)
        x1 <- substr(x, n+1, nchar(x))
        return(paste0(x0, "\n", newline_every_n_chars(x1, n)))
    } else {
        return(x)
    }
}

# Convert to numeric replacing NA with zero
as_numeric_na_zero <- function(x) {
    y <- as.numeric(x)
    y[is.na(y)] <- 0
    return(y)
}

#' @importFrom flextable fontsize autofit dim_pretty width height_all fit_to_width
flextable_fit <- function(ft, width, height,
                          start_size = 11,
                          min_size = 5,
                          max_size = 20,
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



#' @importFrom png readPNG writePNG
#' @importFrom stringr str_c str_remove
crop_png <- function(filename,
                     overwrite = FALSE,
                     output = str_c(str_remove(filename, 'png$'), 'cropped.png'),
                     left = 0, 
                     right = 0) {
    png <- readPNG(filename)
    if (overwrite) { 
        output <- filename 
    }
    writePNG(png[,seq.int(left+1, dim(png)[2] - right),],
             target = output)
}

# read png with png::readPNG to get dimensions
# create external_image with officer
#' @importFrom officer external_img
read_png <- function(filename, dpi = 300) {
    png_dim <- dim(readPNG(filename))[1:2] / dpi
    external_img(filename, width = png_dim[2], height = png_dim[1])
}

remove_child_dirs <- function(dirs) {
    map(dirs, ~ which(str_detect(dirs, str_c('^', .)) &
                          ! str_detect(dirs, str_c('^', ., '$')))) %>%
        unlist() %>% 
        { dirs[-.] }
}

