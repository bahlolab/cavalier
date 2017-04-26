#' Load variants from Bahlo Lab genePanel samplesbed output file
#' 
#' @param genePanel_output_filename genePanel output filename (file ending in [...]_samplesbed.txt)
#' @param sampleID list of IDs and names for samples of interest
# #' @examples
# #' ***TODO***

load_BahloLab_genePanel <- function(genePanel_output_filename, sampleID)
{
    samplesbed <- read.delim(genePanel_output_filename, stringsAsFactors=FALSE)

    # Select and rename relevant columns
    samplesbed_cols <- c("X.Chr", "Start", "End", "Ref", "Alt", "FILTER", "QUAL", "Gene.knownGene", "Func.knownGene", "ExonicFunc.knownGene", "AAChange.knownGene", "snp138", "exac03", "esp6500siv2_all", "X1000g2015aug_all", "SIFT_pred", "Polyphen2_HVAR_pred", "clinvar_20150629")
    vars_cols <- c("chromosome", "start", "end", "reference", "alternate", "FILTER", "QUAL", "gene", "region", "change", "annotation", "dbSNP", "MAF ExAC", "MAF ESP6500", "MAF 1000G", "SIFT", "Polyphen2", "ClinVar")

    vars <- samplesbed[, samplesbed_cols]
    colnames(vars) <- vars_cols

    # Add genotype information for each sampleID
    # (sampleIDs are renamed for ease of understanding and consistency for later trio filtering, etc)
    for (sID in names(sampleID)) {
        sID_name <- sampleID[[sID]]
        if (substring(sID, 1, 1) %in% as.character(0:9)) {
            sID <- paste0("X", sID)
        }
        
        vars[, paste(sID_name, "genotype")] <- samplesbed[, paste0(sID, "_GT")]
        vars[, paste(sID_name, "GT quality")] <- samplesbed[, paste0(sID, "_GQ")]
        vars[, paste(sID_name, "depth (R,A)")] <- samplesbed[, paste0(sID, "_AD")]
    }

    # Replace genes reported twice, separated by a semicolon, to just single gene (eg: ESR1;ESR1 --> ESR1)
    for (ii in 1:nrow(vars)) {
        ii_split_gene <- strsplit(vars[ii, "gene"], ";")[[1]]
        vars[ii, "gene"] <- paste(unique(ii_split_gene), collapse=";")
    }

    # Convert SIFT and Polyphen2 predictions to words
    vars$SIFT[vars$SIFT == "T"] <- "tolerated"
    vars$SIFT[vars$SIFT == "D"] <- "damaging"
    vars$Polyphen2[vars$Polyphen2 == "B"] <- "benign"
    vars$Polyphen2[vars$Polyphen2 == "D"] <- "probably damaging"
    vars$Polyphen2[vars$Polyphen2 == "P"] <- "possibly damaging"    

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
    vars$RVIS <- rvis_exac_percentile(vars$gene)
    vars$Grantham <- as.character(grantham_score(vars$annotation))
    vars$Grantham[is.na(vars$Grantham)] <- ""

    return(vars)
}
