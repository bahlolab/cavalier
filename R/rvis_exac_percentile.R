#' Return RVIS intolerance percentile for given gene symbols
#' 
#' @param gene gene symbols
#' @return RVIS intolerance percentiles
# #' @examples
# #' ***TODO***

rvis_exac_percentile <- function(genes) {
    # RVIS_ExACv2_March2017 table from http://genic-intolerance.org stored as internal data
  miss <- which(!genes %in% RVIS_ExACv2_March2017$gene)
  genes[miss] <- HGNC_alias$symbol[match(genes[miss], HGNC_alias$alias)]
  RVIS_ExACv2_March2017$RVIS_percentile[match(genes, RVIS_ExACv2_March2017$gene)]
}

gevir_percentile <- function(genes) {
  miss <- which(!genes %in% GeVIR$symbol)
  genes[miss] <- HGNC_alias$symbol[match(genes[miss], HGNC_alias$alias)]
  GeVIR$gevir_percentile[match(genes, GeVIR$symbol)]
}

loeuf_percentile <- function(genes) {
  miss <- which(!genes %in% GeVIR$symbol)
  genes[miss] <- HGNC_alias$symbol[match(genes[miss], HGNC_alias$alias)]
  GeVIR$loeuf_percentile[match(genes, GeVIR$symbol)]
}
