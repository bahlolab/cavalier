# Run this script to create internal data tables of public OMIM and RVIS data: R/sysdata.rda

# get approved names and aliases from HGNC
# downloaded from https://www.genenames.org/download/custom/ June 2021
HGNC <- readr::read_tsv('HGNC_2021_06_07.txt.gz', col_types = readr::cols())
HGNC_alias <- dplyr::select(HGNC, symbol = `Approved symbol`)
HGNC_alias$alias <- purrr::map2(HGNC$`Previous symbols`, HGNC$`Alias symbols`,
                                ~ c(stringr::str_split(.x, ', ', simplify = T),
                                    stringr::str_split(.y, ', ', simplify = T)))
HGNC_alias <- tidyr::unnest(HGNC_alias, alias)
# remove ambiguities
HGNC_alias <- dplyr::filter(HGNC_alias, !alias %in% symbol, !is.na(alias))
HGNC_alias <- dplyr::add_count(HGNC_alias, alias)
HGNC_alias <-  dplyr::select(dplyr::filter(HGNC_alias, n == 1),
                             symbol, alias)

hgnc_name_replace <- function(genes) {
  at <- which(genes %in% HGNC_alias$alias)
  replace(genes, at, HGNC_alias$symbol[match(genes[at], HGNC_alias$alias)])
}

# Download RVIS_Unpublished_ExACv2_March2017.txt from http://genic-intolerance.org
RVIS_ExACv2_March2017 <- readr::read_delim("RVIS_Unpublished_ExACv2_March2017.txt.gz", delim="\t")
colnames(RVIS_ExACv2_March2017) <- c("gene", "gene_coverage", "RVIS", "RVIS_percentile", "edge_case_RVIS", "OE_ratio", "OE_ratio_percentile", "alternative_RVIS", "alternative_RVIS_percentile")
RVIS_ExACv2_March2017 <- as.data.frame(RVIS_ExACv2_March2017[, c("gene", "RVIS_percentile")])
RVIS_ExACv2_March2017$gene <- hgnc_name_replace(RVIS_ExACv2_March2017$gene)
rownames(RVIS_ExACv2_March2017) <- RVIS_ExACv2_March2017$gene
# RVIS_ExACv2_March2017$RVIS_percentile <- round(RVIS_ExACv2_March2017$RVIS_percentile, 1)

# Convert between MIM id and gene using public mim2gene table
mim2gene <- read.delim("mim2gene.txt.gz", skip=4, stringsAsFactors=FALSE)
mim2gene[,4] <- hgnc_name_replace(mim2gene[,4])

# note: downloaded June 2021
GeVIR <- readr::read_csv('gevir_gene_rankings.csv.gz', col_types = readr::cols())
GeVIR$symbol <- hgnc_name_replace(GeVIR$symbol)

usethis::use_data(RVIS_ExACv2_March2017, mim2gene, GeVIR, HGNC_alias, internal=TRUE, overwrite=TRUE)

