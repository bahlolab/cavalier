#' Return RVIS intolerance percentile for given gene symbols
#' 
#' @param gene gene symbols
#' @return RVIS intolerance percentiles
# #' @examples
# #' ***TODO***

rvis_exac_percentile <- function(genes) {
  # RVIS_ExACv2_March2017 table from http://genic-intolerance.org stored as internal data
  RVIS_ExACv2_March2017$RVIS_percentile[match(genes, RVIS_ExACv2_March2017$gene)]
}

gevir_percentile <- function(genes) {
  GeVIR$gevir_percentile[match(genes, GeVIR$symbol)]
}

loeuf_percentile <- function(genes) {
  GeVIR$loeuf_percentile[match(genes, GeVIR$symbol)]
}
