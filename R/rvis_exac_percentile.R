#' Return RVIS intolerance percentile for given gene symbols
#' 
#' @param gene gene symbols
#' @return RVIS intolerance percentiles
# #' @examples
# #' ***TODO***

rvis_exac_percentile <- function(genes)
{
    # RVIS_ExACv2_March2017 table from http://genic-intolerance.org stored as internal data
    return(RVIS_ExACv2_March2017[genes, "RVIS_percentile"])
}
