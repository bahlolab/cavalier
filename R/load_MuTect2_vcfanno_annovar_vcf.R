#' Load variants from VCF file (annotations by ANNOVAR, ExAC counts by vcfanno)
#' 
#' @param vcf_filename VEP annotated VCF filename
#' @param sampleID list of IDs and names for samples of interest
# #' @examples
# #' ***TODO***

load_MuTect2_vcfanno_annovar_vcf <- function(vcf_filename, sampleID)
{
    # *** Eventually want the desired list of output columns to be something that can be specified in project settings ***
    vcf_info_fields <- NULL
    vcf_format_fields <- NULL

    vcf <- vcfR::read.vcfR(vcf_filename)
    vcf_tidy <- vcfR::vcfR2tidy(vcf, info_fields=vcf_info_fields, format_fields=vcf_format_fields)

    # Combine and clean VCF data...
    vars <- vcf_tidy$fix[, c("CHROM", "POS", "REF", "ALT", "QUAL", "FILTER", "ExAC_AC", "ExAC_AN", "ExAC_Hom",  "ExAC_FILTER", "gnomAD_exome_AC", "gnomAD_exome_AN", "gnomAD_exome_Hom", "gnomAD_exome_FILTER", "gnomAD_genome_AC", "gnomAD_genome_AN", "gnomAD_genome_Hom", "gnomAD_genome_FILTER", "Gene.refGene", "Func.refGene", "ExonicFunc.refGene", "AAChange.refGene", "avsnp150", "ExAC_ALL", "gnomAD_exome_ALL", "gnomAD_genome_ALL", "SIFT_pred", "Polyphen2_HVAR_pred")]
    colnames(vars) <- c("chromosome", "start", "reference", "alternate", "QUAL", "FILTER", "ExAC count", "ExAC coverage", "ExAC hom count", "ExAC_FILTER", "gnomAD exome count", "gnomAD exome coverage", "gnomAD exome hom count", "gnomAD_exome_FILTER", "gnomAD genome count", "gnomAD genome coverage", "gnomAD genome hom count", "gnomAD_genome_FILTER", "gene", "region", "change", "annotation", "dbSNP", "MAF ExAC", "MAF gnomAD exome", "MAF gnomAD genome", "SIFT", "Polyphen2")
    
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

    vars$change <- gsub("_", " ", vars$change, fixed=TRUE)
    
    # Add missing population frequency columns
    vars$"MAF ESP6500" <- 0
    vars$"MAF 1000G" <- 0
    vars$"MAF ExAC" <- as.numeric(vars$"MAF ExAC")
    vars$"MAF ExAC"[is.na(vars$"MAF ExAC")] <- 0
    
    vars$"gnomAD exome count" <- as.numeric(vars$"gnomAD exome count")
    vars$"gnomAD genome count" <- as.numeric(vars$"gnomAD genome count")
    vars$"gnomAD exome count"[is.na(vars$"gnomAD exome count")] <- 0
    vars$"gnomAD genome count"[is.na(vars$"gnomAD genome count")] <- 0
    vars$"gnomAD count" <- vars$"gnomAD exome count" + vars$"gnomAD genome count"
    
    vars$"gnomAD exome hom count" <- as.numeric(vars$"gnomAD exome hom count")
    vars$"gnomAD genome hom count" <- as.numeric(vars$"gnomAD genome hom count")
    vars$"gnomAD exome hom count"[is.na(vars$"gnomAD exome hom count")] <- 0
    vars$"gnomAD genome hom count"[is.na(vars$"gnomAD genome hom count")] <- 0
    vars$"gnomAD hom count" <- vars$"gnomAD exome hom count" + vars$"gnomAD genome hom count"
    
    # Convert SIFT and Polyphen2 predictions to words
    vars$SIFT[vars$SIFT == "T"] <- "tolerated"
    vars$SIFT[vars$SIFT == "D"] <- "damaging"
    vars$Polyphen2[vars$Polyphen2 == "B"] <- "benign"
    vars$Polyphen2[vars$Polyphen2 == "D"] <- "probably damaging"
    vars$Polyphen2[vars$Polyphen2 == "P"] <- "possibly damaging"    
    
    vars$RVIS <- rvis_exac_percentile(vars$gene)
    vars$Grantham <- as.character(grantham_score(vars$annotation))
    vars$Grantham[is.na(vars$Grantham)] <- ""

    return(as.data.frame(vars))
}

