#' Load variants from Bahlo Lab MPSPipeline merged_collab output file
#' 
#' @param MPSPipeline_output_filename MPSPipeline output filename (file ending in [...]_merged_collab.txt)
#' @param sampleID list of IDs and names for samples of interest
# #' @examples
# #' ***TODO***

load_BahloLab_MPSPipeline <- function(MPSPipeline_output_filename, sampleID)

{
    merged_collab <- read.delim(MPSPipeline_output_filename, stringsAsFactors=FALSE)

    merged_collab_cols <- c("CHROM", "START", "END", "REF", "ALT", "Filter", "QUAL", "gene", "region", "change", "annotation", "dbSNP135", "ExAC", "esp6500_all", "aaf.1KG", "SIFT.Pred", "pph2_pred")
    vars_cols <- c("chromosome", "start", "end", "reference", "alternate", "FILTER", "QUAL", "gene", "region", "change", "annotation", "dbSNP", "MAF ExAC", "MAF ESP6500", "MAF 1000G", "SIFT", "Polyphen2")

    vars <- merged_collab[, merged_collab_cols]
    colnames(vars) <- vars_cols

    vars$"inhouse control AC" <- as.numeric(gsub(".", "0", merged_collab$control.AC, fixed=TRUE))

    # Add genotype information for each sampleID
    # (sampleIDs are renamed for ease of understanding and consistency for later trio filtering, etc)
    for (sID in names(sampleID)) {
        sID_name <- sampleID[[sID]]
        if (substring(sID, 1, 1) %in% as.character(0:9)) {
            sID <- paste0("X", sID)
        }
        vars[, paste(sID_name, "genotype")] <- merged_collab[, paste0(sID, "_GT")]
        vars[, paste(sID_name, "GT quality")] <- merged_collab[, paste0(sID, "_GQ")]
        vars[, paste(sID_name, "depth (R,A)")] <- paste0("(", merged_collab[, paste0(sID, "_DPR")], ",", merged_collab[, paste0(sID, "_DPA")], ")")

        # Convert genotype format from "1", "2", "2_1/2" to "0/1", "1/1", "1/2" respectively
        # (I prefer this format aesthetically and it is consistent with genePanel results below)
        sID_GT <- vars[, paste(sID_name, "genotype")]
        sID_GT <- gsub("1_", "", gsub("2_", "", sID_GT))
        sID_GT <- ifelse(sID_GT == "0", "0/0", sID_GT)
        sID_GT <- ifelse(sID_GT == "1", "0/1", sID_GT)
        sID_GT <- ifelse(sID_GT == "2", "1/1", sID_GT)
        sID_GT <- ifelse(sID_GT == ".", NA, sID_GT)
        vars[, paste(sID_name, "genotype")] <- sID_GT
    }

    # Convert SIFT format
    vars$SIFT <- tolower(vars$SIFT)

    # Splicing vars have GENE(annotation) stored in gene column, so parse it and separate to gene and annotation columns
    parse_gene <- sapply(vars$gene, function(x){strsplit(x, "(", fixed=TRUE)[[1]][1]})
    parse_ann <- sapply(vars$gene, function(x){gsub(")", "", strsplit(x, "(", fixed=TRUE)[[1]][2], fixed=TRUE)})
    vars$gene <- parse_gene
    vars$annotation <- ifelse(vars$annotation == "." & !is.na(parse_ann), paste(parse_gene, parse_ann, sep=":"), vars$annotation)

    # Replace genes reported twice, separated by a semicolon, to just single gene (eg: ESR1;ESR1 --> ESR1)
    for (ii in 1:nrow(vars)) {
        ii_split_gene <- strsplit(vars[ii, "gene"], ";")[[1]]
        vars[ii, "gene"] <- paste(unique(ii_split_gene), collapse=";")
    }

    # Convert MAF columns to numeric values with "." equal to 0
    vars$`MAF ExAC`[vars$`MAF ExAC` == "."] <- "0"
    vars$`MAF ExAC` <- as.numeric(vars$`MAF ExAC`)
    vars$`MAF ESP6500`[vars$`MAF ESP6500` == "."] <- "0"
    vars$`MAF ESP6500` <- as.numeric(vars$`MAF ESP6500`)
    vars$`MAF 1000G`[vars$`MAF 1000G` == "."] <- "0"
    vars$`MAF 1000G` <- as.numeric(vars$`MAF 1000G`)

    # *** Temporary hack to get approximate ExAC count (as compared to ExAC MAF) ***
    vars[, "ExAC count"] <- round(120000 * vars[, "MAF ExAC"])

    # Add RVIS ExAC score (http://genic-intolerance.org) and Grantham scores
    vars$RVIS <- RVIS_ExAC_percentile(vars$gene)
    vars$Grantham <- as.character(Grantham_score(vars$annotation))
    vars$Grantham[is.na(vars$Grantham)] <- ""
    
    return(vars)
}

