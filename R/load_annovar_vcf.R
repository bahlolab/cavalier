#' Load variants from VCF file
#' 
#' @param vcf_filename ANNOVAR annotated VCF filename
#' @param sampleID list of IDs and names for samples of interest
#' @param rename_columns column names in VCF and names to rename them to
#' @param numeric_columns columns to convert to numeric values
# #' @examples
# #' ***TODO***

load_annovar_vcf <- function(vcf_filename, sampleID, rename_columns=NULL, numeric_columns=NULL)
{
    vcf_info_fields <- NULL
    vcf_format_fields <- NULL

    vcf <- vcfR::read.vcfR(vcf_filename)
    vcf_tidy <- vcfR::vcfR2tidy(vcf, info_fields=vcf_info_fields, format_fields=vcf_format_fields)

    vars <- as.data.frame(vcf_tidy$fix)
    
    for (sID in names(sampleID)) {
        vcf_tidy_gt_sID <- vcf_tidy$gt[vcf_tidy$gt$Indiv == sID, ]

        sID_name <- sampleID[[sID]]
        vars[, paste(sID_name, "genotype")] <- vcf_tidy_gt_sID[, "gt_GT"]
        vars[, paste(sID_name, "GT quality")] <- vcf_tidy_gt_sID[, "gt_GQ"]
        vars[, paste(sID_name, "depth (R,A)")] <- vcf_tidy_gt_sID[, "gt_AD"]
    }
    
    # Force all chromosomes to start with "chr", eg "chr1" rather than "1"
    vars[, 2] <- ifelse(startsWith(vars[, 2], "chr"), vars[, 2], paste0("chr", vars[, 2]))
    
    # Default column renaming
    if (is.null(rename_columns)) {
        rename_columns <- c("CHROM"="chromosome",
                            "POS"="position",
                            "REF"="reference",
                            "ALT"="alternate",
                            "Gene.refGene"="gene",
                            "Func.refGene"="region",
                            "ExonicFunc.refGene"="change",
                            "AAChange.refGene"="annotation",
                            "SIFT_pred"="SIFT",
                            "Polyphen2_HVAR_pred"="Polyphen2",
                            "CADD_phred"="CADD",
                            "ExAC_AC"="ExAC count",
                            "ExAC_AN"="ExAC coverage",
                            "ExAC_Hom"="ExAC hom count",
                            "gnomAD_exome_AC"="gnomAD exome count",
                            "gnomAD_exome_AN"="gnomAD exome coverage",
                            "gnomAD_exome_Hom"="gnomAD exome hom count",
                            "gnomAD_genome_AC"="gnomAD genome count",
                            "gnomAD_genome_AN"="gnomAD genome coverage",
                            "gnomAD_genome_Hom"="gnomAD genome hom count")
    }
    
    for (cc in colnames(vars)[colnames(vars) %in% names(rename_columns)]) {
        colnames(vars)[which(colnames(vars) == cc)] <- rename_columns[cc]
    }
    
    # Fix ANNOVAR character format
    fix_format_columns <- c("gene", "Gene.refGene", "GeneDetail.refGene", "region", "ExonicFunc.refGene", "annotation", "AAChange.refGene", "genomicSuperDups")
    for (cc in fix_format_columns) {
        if (cc %in% colnames(vars)) {
            vars[, cc] <- gsub("\\x3b", ";", gsub("\\x3d", "=", vars[, cc], fixed=TRUE), fixed=TRUE)
        }
    }
    
    # Split annotation and select longest only the longest transcript
    # *** FUNCTION NEEDS WORK ***
    get_annotation_longest_transcript <- function(annotations) {
        annot_split <- unlist(strsplit(annotations, ","))
        transcript <- unlist(lapply(annot_split, function(x) {y <- strsplit(x, ":")[[1]]; y[2]}))
        codon <- unlist(lapply(annot_split, function(x) {y <- strsplit(x, ":")[[1]]; y[startsWith(y, "c.")]}))
        prot <- unlist(lapply(annot_split, function(x) {y <- strsplit(x, ":")[[1]]; y[startsWith(y, "p.")]}))
        if (length(codon) == 0) {
            return(annotations)
        }
        bp_num <- sapply(codon, function(x){y <- strsplit(x, "")[[1]]; as.integer(paste(y[y %in% as.character(0:9)], collapse=""))})
        longest <- which(bp_num == max(bp_num))[1]
        # return(paste(transcript[longest], codon[longest], prot[longest], sep=":"))
        return(paste(codon[longest], prot[longest], sep=":"))
    }
    
    if ("annotation" %in% colnames(vars)) {
        vars$annotation <- sapply(vars$annotation, get_annotation_longest_transcript)
        vars$Grantham <- as.character(grantham_score(vars$annotation))
        vars$Grantham[is.na(vars$Grantham)] <- ""
    }
    
    if ("gene" %in% colnames(vars)) {
        vars$RVIS <- rvis_exac_percentile(vars$gene)
    }
    
    # Modify region and change for splicing variants (so synonymous SNV;splicing variants are not removed with other synonymous SNV)
    # *** TODO: clean this up and test
    if (all(c("change", "region") %in% colnames(vars))) {
        vars$change[vars$region == "splicing"] <- "splicing"
        vars$change[vars$region == "exonic;splicing"] <- paste0(vars$change[vars$region == "exonic;splicing"], ";splicing")
    }
    # Collapse duplicate gene names for variants in exonic and splice region
    if (all(c("gene", "region") %in% colnames(vars))) {
        vars$gene[vars$region == "exonic;splicing"] <- sapply(strsplit(vars$gene[vars$region == "exonic;splicing"], "\\", fixed=TRUE), function(x){x[1]})
    }
    
    if (all(c("gnomAD exome count", "gnomAD genome count") %in% colnames(vars))) {
        vars$"gnomAD count" <- as_numeric_na_zero(vars$"gnomAD exome count") + as_numeric_na_zero(vars$"gnomAD genome count")
    }
    if (all(c("gnomAD exome hom count", "gnomAD genome hom count") %in% colnames(vars))) {
        vars$"gnomAD hom count" <- as_numeric_na_zero(vars$"gnomAD exome hom count") + as_numeric_na_zero(vars$"gnomAD genome hom count")
    }
    
    return(vars)
}

