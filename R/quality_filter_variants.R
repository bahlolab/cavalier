#' Filter variants based on quality scores (use GATK best practices as default)
#' 
#' @param vars variants data.frame
#' @return remaining variants data.frame after applying quality filters
# #' @examples
# #' ***TODO***

# *** Use default values copied from Broad Institue GATK Best Practices but remove this if using VQSR ***
# *** In future want to be able to specify quality filter settings in project settings ***
# 
# https://software.broadinstitute.org/gatk/guide/article?id=3225
# 
# For SNPs:
# QD < 2.0
# MQ < 40.0
# FS > 60.0
# SOR > 3.0
# MQRankSum < -12.5
# ReadPosRankSum < -8.0
# 
# For indels:
# QD < 2.0
# ReadPosRankSum < -20.0
# InbreedingCoeff < -0.8
# FS > 200.0
# SOR > 10.0

quality_filter_variants <- function(vars)
{
    isSNV <- (nchar(vars$reference) == 1) & (nchar(vars$alternate) == 1)

    if ("QD" %in% colnames(vars)) {
        QD_remove <- vars$QD < 2
    } else {
        QD_remove <- rep(FALSE, nrow(vars))
    }
    if ("MQ" %in% colnames(vars)) {
        MQ_remove <- ifelse(isSNV, vars$MQ < 40, FALSE)
    } else {
        MQ_remove <- rep(FALSE, nrow(vars))
    }
    if ("FS" %in% colnames(vars)) {
        FS_remove <- ifelse(isSNV, vars$FS > 60, vars$FS > 200)
    } else {
        FS_remove <- rep(FALSE, nrow(vars))
    }
    if ("SOR" %in% colnames(vars)) {
        SOR_remove <- ifelse(isSNV, vars$SOR > 3, vars$SOR > 10)
    } else {
        SOR_remove <- rep(FALSE, nrow(vars))
    }
    if ("MQRankSum" %in% colnames(vars)) {
        MQRankSum_remove <- ifelse(isSNV, vars$MQRankSum < -12.5, FALSE)
    } else {
        MQRankSum_remove <- rep(FALSE, nrow(vars))
    }
    if ("ReadPosRankSum" %in% colnames(vars)) {
        ReadPosRankSum_remove <- ifelse(isSNV, vars$ReadPosRankSum < -8, vars$ReadPosRankSum < -20)
    } else {
        ReadPosRankSum_remove <- rep(FALSE, nrow(vars))
    }
    if ("InbreedingCoeff" %in% colnames(vars)) {
        InbreedingCoeff_remove <- ifelse(isSNV, FALSE, vars$InbreedingCoeff < -0.8)
    } else {
        InbreedingCoeff_remove <- rep(FALSE, nrow(vars))
    }
    
    quality_remove <- QD_remove | MQ_remove | FS_remove | SOR_remove | MQRankSum_remove | ReadPosRankSum_remove | InbreedingCoeff_remove
    quality_remove[is.na(quality_remove)] <- FALSE

    filter_vars <- vars[!quality_remove, ]

    return(filter_vars)
}

