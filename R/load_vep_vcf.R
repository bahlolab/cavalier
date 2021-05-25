#' Load variants from VCF file (ANN field from VEP, ExAC counts by vcfanno)
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
    
    # Get ANN format column names
    ANN_columns <- strsplit(strsplit(vcf_tidy$meta$Description[vcf_tidy$meta$ID == "ANN"], "Format: ")[[1]][2], "|", fixed=TRUE)[[1]]
    
    ANN_split1 <- strsplit(vcf_tidy$fix$ANN, ",")
    
    get_canonical_ANN <- function(ANN) {
        ANN_df <- as.data.frame(do.call(rbind, strsplit(ANN, "|", fixed=TRUE)), stringsAsFactors=FALSE)
        colnames(ANN_df) <- ANN_columns[1:ncol(ANN_df)]
        ANN_df <- ANN_df[ANN_df$CANONICAL == "YES", ]
        apply(ANN_df, 2, function(x){paste(unique(x), collapse=",")})
    }
    
    vcf_tidy[["ANN"]] <- as.data.frame(do.call(rbind, lapply(ANN_split1, get_canonical_ANN)), stringsAsFactors=FALSE)
    colnames(vcf_tidy$ANN) <- ANN_columns
    rownames(vcf_tidy$ANN) <- rownames(vcf_tidy$fix)
    
    # Combine and clean VCF data
    vars <- as.data.frame(vcf_tidy$fix[, c("CHROM", "POS", "REF", "ALT", "QUAL", "FILTER", "DP", "QD", "MQ", "FS", "SOR", "MQRankSum", "ReadPosRankSum", "InbreedingCoeff")], stringsAsFactors=FALSE)
    vars <- cbind(vars, vcf_tidy$ANN[, c("Consequence", "SYMBOL", "Existing_variation", "SIFT", "PolyPhen", "AF", "gnomAD_AF")])

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
    AA_split <- strsplit(vcf_tidy$ANN$Amino_acids, "/", fixed=TRUE)
    AA1 <- sapply(AA_split, function(x){x[1]})
    AA2 <- sapply(AA_split, function(x){x[2]})
    AA1[is.na(AA1)] <- ""
    AA2[is.na(AA2)] <- ""
    vars$annotation <- paste0("p.", AA1, vcf_tidy$ANN$Protein_position, AA2)
    vars$annotation <- gsub(",", "", vars$annotation, fixed=TRUE)
    vars$annotation[vars$annotation == "p."] <- ""
    
    vars$dbSNP <- sapply(strsplit(vars$dbSNP, "&"), function(x){paste(x[startsWith(x, "rs")], collapse=";")})

    # *** TODO: add more population frequency information and exact ExAC count ***
    vars$"MAF 1000G" <- as.numeric(vars$"MAF 1000G")
    vars$"MAF 1000G"[is.na(vars$"MAF 1000G")] <- 0
    vars$"MAF ExAC" <- as.numeric(vars$"MAF ExAC")
    vars$"MAF ExAC"[is.na(vars$"MAF ExAC")] <- 0
    

    # sift_regex <- "([a-z]*)\\(([0-9\\.]*)\\)"
    # vars$SIFT <- sub(sift_regex, "\\1", vars$SIFT)
    # vars$SIFT <- gsub("_", " ", gsub("_low_confidence", " *warning! low confidence.", vars$SIFT))
    # polyphen_regex <- "([a-z_]*)\\(([0-9\\.]*)\\)"
    # vars$Polyphen2 <- sub(polyphen_regex, "\\1", vars$Polyphen2)
    # vars$Polyphen2 <- gsub("_", " ", vars$Polyphen2)

    vars$RVIS <- rvis_exac_percentile(vars$gene)
    vars$Grantham <- grantham_score(vars$annotation)

    return(vars)
}

