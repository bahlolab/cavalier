#' Load variants from VCF file (CSQ field from VEP, ExAC counts by vcfanno)
#' 
#' @param vcf_filename VEP annotated VCF filename
#' @param sampleID list of IDs and names for samples of interest
# #' @examples
# #' ***TODO***

load_vep_vcf <- function(vcf_filename, sampleID)
{
    # *** Eventually want the desired list of output columns to be something that can be specified in project settings ***
    vcf_info_fields <- NULL
    vcf_format_fields <- NULL

    vcf <- vcfR::read.vcfR(vcf_filename)
    vcf_tidy <- vcfR::vcfR2tidy(vcf, info_fields=vcf_info_fields, format_fields=vcf_format_fields)

    # *** Currently take only first CSQ but need to handle cases where multiple are present more robustly ***
    CSQ_columns <- strsplit(sub("\">", "", strsplit(vcf@meta[startsWith(vcf@meta, "##INFO=<ID=CSQ")], "Format: ")[[1]][2]), "|", fixed=TRUE)[[1]]
    vcf_tidy[["CSQ"]] <- as.data.frame(do.call(rbind, lapply(strsplit(sapply(strsplit(vcf_tidy$fix$CSQ, ",", fixed=TRUE), function(x)x[1]), "|", fixed=TRUE), function(x){x[1:length(CSQ_columns)]})), stringsAsFactors=FALSE)
    colnames(vcf_tidy$CSQ) <- CSQ_columns
    rownames(vcf_tidy$CSQ) <- rownames(vcf_tidy$fix)
    

    # Combine and clean VCF data...
    vars <- vcf_tidy$fix[, c("CHROM", "POS", "REF", "ALT", "QUAL", "FILTER", "DP", "QD", "MQ", "FS", "SOR", "MQRankSum", "ReadPosRankSum", "InbreedingCoeff")]
    vars <- cbind(vars, vcf_tidy$CSQ[, c("Consequence", "SYMBOL", "Existing_variation", "SIFT", "PolyPhen", "GMAF", "ExAC_Adj_MAF")])

    colnames(vars) <- c("chromosome", "start", "reference", "alternate", "QUAL", "FILTER", "DP", "QD", "MQ", "FS", "SOR", "MQRankSum", "ReadPosRankSum", "InbreedingCoeff", "change", "gene", "dbSNP", "SIFT", "Polyphen2", "MAF 1000G", "MAF ExAC")

    for (sID in names(sampleID)) {
        vcf_tidy_gt_sID <- vcf_tidy$gt[vcf_tidy$gt$Indiv == sID, ]
        # *** BETTER TESTING SOON ***
        stopifnot(all(vcf_tidy_gt_sID$POS == vars$start))

        sID_name <- sampleID[[sID]]
        vars[, paste(sID_name, "genotype")] <- vcf_tidy_gt_sID[, "gt_GT"]
        vars[, paste(sID_name, "GT quality")] <- vcf_tidy_gt_sID[, "gt_GQ"]
        vars[, paste(sID_name, "depth (R,A)")] <- vcf_tidy_gt_sID[, "gt_AD"]
    }

    vars$chromosome <- ifelse(startsWith(vars$chromosome, "chr"), vars$chromosome, paste0("chr", vars$chromosome))

    # Construct end position from start position and difference between length of ref and alt alleles
    length_ref <- sapply(vars$reference, nchar)
    length_alt <- sapply(vars$alternate, nchar)
    # Start and end positions are the same for insertions or SNV but different for deletions
    vars$end <- ifelse(length_ref <= length_alt, vars$start, vars$start + length_ref - length_alt - 1)

    vars$change <- gsub("&", ";", vars$change, fixed=TRUE)
    
    # Annotation  *** TODO: add
    vars$annotation <- ""
    
    vars$dbSNP <- sapply(strsplit(vars$dbSNP, "&"), function(x){paste(x[startsWith(x, "rs")], collapse=";")})

    # *** TODO: add more population frequency information and exact ExAC count ***
    vars$"MAF 1000G" <- as.numeric(sapply(strsplit(sapply(strsplit(vars$"MAF 1000G", "&"), function(x)x[1]), ":"), function(x)x[2]))
    vars$"MAF ExAC" <- as.numeric(sapply(strsplit(sapply(strsplit(vars$"MAF ExAC", "&"), function(x)x[1]), ":"), function(x)x[2]))
    vars$"MAF ESP6500" <- as.numeric(".")
    vars$"ExAC count" <- 120000 * as.numeric(vars$"MAF ExAC")
    vars$"ExAC count"[is.na(vars$"ExAC count")] <- 0
    vars$"MAF ExAC"[is.na(vars$"MAF ExAC")] <- 0
    vars$"MAF 1000G"[is.na(vars$"MAF 1000G")] <- 0
    vars$"MAF ESP6500"[is.na(vars$"MAF ESP6500")] <- 0

    sift_regex <- "([a-z]*)\\(([0-9\\.]*)\\)"
    vars$SIFT <- sub(sift_regex, "\\1", vars$SIFT)
    vars$SIFT <- gsub("_", " ", gsub("_low_confidence", " *warning! low confidence.", vars$SIFT))
    polyphen_regex <- "([a-z_]*)\\(([0-9\\.]*)\\)"
    vars$Polyphen2 <- sub(polyphen_regex, "\\1", vars$Polyphen2)
    vars$Polyphen2 <- gsub("_", " ", vars$Polyphen2)

    vars$RVIS <- rvis_exac_percentile(vars$gene)
    vars$Grantham <- ""  # *** NEED TO FIX GRANTHAM FOR DIFFERENT ANNOTATION

    return(vars)
}

