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

# replace gene symbols alias with HGNC approved symbol
#' @export
hgnc_name_replace <- function(genes) {
    at <- which(genes %in% HGNC_alias$alias)
    replace(genes, at, HGNC_alias$symbol[match(genes[at], HGNC_alias$alias)])
}

