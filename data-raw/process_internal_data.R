# Run this script to create internal data tables of public OMIM and RVIS data: R/sysdata.rda

# Download RVIS_Unpublished_ExACv2_March2017.txt from http://genic-intolerance.org
RVIS_ExACv2_March2017 <- readr::read_delim("RVIS_Unpublished_ExACv2_March2017.txt", delim="\t")
colnames(RVIS_ExACv2_March2017) <- c("gene", "gene_coverage", "RVIS", "RVIS_percentile", "edge_case_RVIS", "OE_ratio", "OE_ratio_percentile", "alternative_RVIS", "alternative_RVIS_percentile")
RVIS_ExACv2_March2017 <- as.data.frame(RVIS_ExACv2_March2017[, c("gene", "RVIS_percentile")])
rownames(RVIS_ExACv2_March2017) <- RVIS_ExACv2_March2017$gene
RVIS_ExACv2_March2017$RVIS_percentile <- round(RVIS_ExACv2_March2017$RVIS_percentile, 1)

# Convert between MIM id and gene using public mim2gene table
mim2gene <- read.delim("mim2gene.txt", skip=4, stringsAsFactors=FALSE)

# Store as internal dataset
devtools::use_data(RVIS_ExACv2_March2017, mim2gene, internal=TRUE, overwrite=TRUE)

