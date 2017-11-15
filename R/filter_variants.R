#' Filter list of variants based on inheritance model, variant type and frequency
#' 
#' @param vars variants data.frame
#' @param inheritance_MAF list of inheritance models and MAF thresholds to filter
#' @param sampleID list of IDs and names for samples of interest
#' @param gene_include genes to include in filtering
#' @param gene_exclude genes to exclude from filtering (if gene_include not specified)
#' @param region_include region to include in filtering
#' @param region_exclude region to exclude from filtering (if region_include not specified)
#' @param change_include change to include in filtering
#' @param change_exclude change to exclude from filtering (if change_include not specified)
#' @param min_depth minimum depth to require in all individuals to keep variant (default: 0)
#' @return data frame of candidate variants for each inheritance model
# #' @examples
# #' ***TODO***

filter_variants <- function(vars, inheritance_MAF, sampleID,
                            gene_include=NULL,  gene_exclude=NULL,
                            region_include=NULL, region_exclude=NULL,
                            change_include=NULL, change_exclude=NULL,
                            min_depth=0)
{
    names(inheritance_MAF) <- gsub("compound heterozygous", "comp het", names(inheritance_MAF))

    # Filter based on genes
    if (!is.null(gene_include)) {
        gene_filter <- toupper(vars$gene) %in% toupper(gene_include)
    } else if (!is.null(gene_exclude)) {
        gene_filter <- !(toupper(vars$gene) %in% toupper(gene_exclude))
    } else {
        gene_filter <- rep(TRUE, nrow(vars))
    }
    # Filter based on region
    if (!is.null(region_include)) {
        region_filter <- vars$region %in% region_include
    } else if (!is.null(region_exclude)) {
        region_filter <- !(vars$region %in% region_exclude)
    } else {
        region_filter <- rep(TRUE, nrow(vars))
    }
    # Filter based on change
    if (!is.null(change_include)) {
        change_filter <- vars$change %in% change_include
    } else if (!is.null(change_exclude)) {
        change_filter <- !(vars$change %in% change_exclude)
    } else {
        change_filter <- rep(TRUE, nrow(vars))
    }

    # Filter for minimum depth (require >= min_depth in all samples)
    depth_filter <- rep(TRUE, nrow(vars))
    if (min_depth > 0) {
        depth_cols <- colnames(vars)[endsWith(colnames(vars), "depth (R,A)")]
        for (dc in depth_cols) {
            dc_depth_col <- gsub("(", "", gsub(")", "", vars[, dc], fixed=TRUE), fixed=TRUE)
            dc_split <- strsplit(dc_depth_col, ",")
            dc_ref <- as.integer(sapply(dc_split, function(x)x[1]))
            dc_alt <- as.integer(sapply(dc_split, function(x)x[2]))
            dc_depth <- dc_ref + dc_alt

            depth_filter <- depth_filter & (dc_depth >= min_depth)
        }
    }

    # Filter
    filter_vars <- vars[gene_filter & region_filter & change_filter & depth_filter, ]

    # Filter inheritance models
    candvars <- NULL
    for (inh_model in names(inheritance_MAF)) {

        new_candvars <- NULL

        # Filter for MAF
        MAF_thresh <- inheritance_MAF[[inh_model]]
        MAF_ExAC_filter <- filter_vars$`MAF ExAC` <= MAF_thresh
        MAF_ESP6500_filter <- filter_vars$`MAF ESP6500` <= MAF_thresh
        MAF_1000G_filter <- filter_vars$`MAF 1000G` <= MAF_thresh
        MAF_filter <- MAF_ExAC_filter & MAF_ESP6500_filter & MAF_1000G_filter
        MAF_filter_vars <- filter_vars[MAF_filter, ]

        # For individual inheritance model use "proband" sample (works even if multiple samples present)
        # or if only one sample listed then  use it regardless of name
        if (startsWith(inh_model, "individual")) {
            if ("proband" %in% sampleID) {
                individual_name <- "proband"
            } else {
                if (length(sampleID) != 1) {
                    stop("Error: invalid sampleID for individual analysis")
                }
                individual_name <- sampleID[[1]]
            }
            individual_GT <- paste(individual_name, "genotype")
        } else if (startsWith(inh_model, "trio") | endsWith(inh_model, "isodisomy")) {
            if (!all(c("proband", "father", "mother") %in% sampleID)) {
                stop("Error: invalid sampleID for trio analysis, names should be proband, father, mother")
            }
            proband_GT <- MAF_filter_vars$`proband genotype`
            father_GT <- MAF_filter_vars$`father genotype`
            mother_GT <- MAF_filter_vars$`mother genotype`
        }

        # *** Need to check genotypes of X and Y chromosome genes more carefully ***

        if (inh_model == "none") {
            new_candvars <- MAF_filter_vars
        }
        else if (inh_model == "individual dominant") {
            new_candvars <- MAF_filter_vars[MAF_filter_vars[, individual_GT] %in% c("0/1"), ]
        }
        else if (inh_model == "individual recessive") {
            new_candvars <- MAF_filter_vars[MAF_filter_vars[, individual_GT] %in% c("1/1"), ]
        }
        else if (inh_model == "individual comp het") {
            onehit <- MAF_filter_vars[MAF_filter_vars[, individual_GT] %in% c("0/1"), ]
            twohit_genes <- names(table(onehit$gene))[table(onehit$gene) >= 2]
            new_candvars <- onehit[onehit$gene %in% twohit_genes, ]
        }
        else if (inh_model == "shared dominant") {
            dominant_in_all_sampleIDs <- rep(TRUE, nrow(MAF_filter_vars))
            for (sID in names(sampleID)) {
                sID_GT <- MAF_filter_vars[, paste(sampleID[[sID]], "genotype")] %in% c("0/1", NA)
                dominant_in_all_sampleIDs <- dominant_in_all_sampleIDs & sID_GT
            }
            new_candvars <- MAF_filter_vars[dominant_in_all_sampleIDs, ]
        }
        else if (inh_model == "trio de novo") {
            new_candvars <- MAF_filter_vars[proband_GT %in% c("0/1") & father_GT %in% c("0/0", NA) & mother_GT %in% c("0/0", NA), ]
        }
        else if (inh_model == "trio recessive") {
            new_candvars <- MAF_filter_vars[proband_GT %in% c("1/1") & father_GT %in% c("0/1", NA) & mother_GT %in% c("0/1", NA), ]
        }
        else if (inh_model == "trio comp het") {

            comphet_father <- MAF_filter_vars[proband_GT %in% c("0/1") & father_GT %in% c("0/1", NA) & mother_GT %in% c("0/0", NA), ]
            comphet_mother <- MAF_filter_vars[proband_GT %in% c("0/1") & mother_GT %in% c("0/1", NA) & father_GT %in% c("0/0", NA), ]
            
            twohits_father <- comphet_father[comphet_father$gene %in% comphet_mother$gene, ]
            twohits_mother <- comphet_mother[comphet_mother$gene %in% comphet_father$gene, ]

            new_candvars <- merge(twohits_father, twohits_mother, all.x=TRUE, all.y=TRUE)
        }
        else if (inh_model == "paternal isodisomy") {
            
            # *** TODO: ADD FILTER TO SPECIFIC REGION ***

            # Potential paternal isodisomy variants
            new_candvars <- MAF_filter_vars[father_GT %in% "0/1" & proband_GT %in% "0/1" & mother_GT %in% c("0/0", "0/1"), ]

            # Filter candidates with alternate > 50% of counts (i.e., they got copy of alt from paternal side)
            proband_depth <- gsub("(", "", gsub(")", "", new_candvars$`proband depth (R,A)`, fixed=TRUE), fixed=TRUE)
            proband_split <- strsplit(proband_depth, ",")
            proband_ref_counts <- as.integer(sapply(proband_split, function(x)x[1]))
            proband_alt_counts <- as.integer(sapply(proband_split, function(x)x[2]))
            proband_alt_allele_freq <- proband_alt_counts / (proband_ref_counts + proband_alt_counts)

            new_candvars <- new_candvars[proband_alt_allele_freq > 0.5, ]
        }
        else if (inh_model %in% c("trio x-linked dominant", "paternal heterodisomy", "maternal isodisomy", "maternal heterodisomy")) {
            stop(paste("Error:", inh_model, "inheritance model not yet implemented"))
        }
        else {
            stop(paste("Error: invalid inheritance model", inh_model))
        }

        if (nrow(new_candvars) > 0) {
            new_candvars[, "inheritance model"] <- inh_model
            candvars <- rbind(candvars, new_candvars)
        }
    }

    return(candvars)
}