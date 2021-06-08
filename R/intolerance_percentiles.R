#' Return RVIS intolerance percentile for given gene symbols
#' 
#' @param gene gene symbols
#' @return RVIS intolerance percentiles
# #' @examples
# #' ***TODO***

rvis_exac_percentile <- function(genes, round = 1) {
  # RVIS_ExACv2_March2017 table from http://genic-intolerance.org stored as internal data
  p <- RVIS_ExACv2_March2017$RVIS_percentile[match(genes, RVIS_ExACv2_March2017$gene)]
  round(p, round)
}

gevir_percentile <- function(genes, round = 1) {
  p <- GeVIR$gevir_percentile[match(genes, GeVIR$gene_name)]
  round(p, round)
}

loeuf_percentile <- function(genes, round = 1) {
  p <- GeVIR$loeuf_percentile[match(genes, GeVIR$gene_name)]
  round(p, round)
}
