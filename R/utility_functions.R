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

#' @importFrom rlang is_double is_integer
is_number <- function(x) { (is_double(x) | is_integer(x)) & !any(is.na(x)) }


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

#' @importFrom purrr walk map
#' @export
clear_cache <- function(mem = TRUE, disk = FALSE)
{
    if (disk) {
        cache_dir <- get_cavalier_opt('cache_dir')
        if (dir.exists(cache_dir)) {
            file.remove(list.files(cache_dir, full.names = TRUE))
        }
    }
    if (mem) {
        rm(list = ls(envir = cavalier_cache), envir = cavalier_cache)
    }
}

# execute a function and save to disk as filename unless filename exists, then load from file
# useful to cache downloaded files
cache <- function(fun, filename) 
{
    cache_dir <- get_cavalier_opt('cache_dir')
    cache_file <- file.path(cache_dir, filename)
    if (!str_ends(cache_file, '.rds')) {
        cache_file <- str_c(cache_file, '.rds')
    }
    
    if (file.exists(cache_file)) {
        readRDS(cache_file)
    } else {
        res <- fun()
        tmp_fn <- tempfile(pattern = basename(cache_file) %>% str_c('.'),
                           tmpdir = cache_dir)
        saveRDS(res, tmp_fn)
        file.rename(tmp_fn, cache_file)
        res
    }
}
# predicates for checking arguments

is_null_or_file <- function(x) { is.null(x) | (is_scalar_character(x) && file.exists(x)) }

is_null_or_files <- function(x, named = FALSE) { 
    is.null(x) | (is_character(x) && all(file.exists(x))) && (!named | is_named(x))
}

assert_create_dir <- function(dirname, recursive = TRUE) {
    if (!dir.exists(dirname)) {
        assert_that(dir.create(dirname, recursive = recursive))
    }
}

