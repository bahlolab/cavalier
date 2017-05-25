#' Create HTML table of candidate variants
#' 
#' @param candidates candidate variants data.frame
#' @param output_dir cavalier output directory
#' @param sampleID list of IDs and names for samples of interest
#' @param GTEx_median_rpkm location of GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct.gz file downloaded from GTEx Portal (see https://gtexportal.org/home/datasets)
#' @param hide_missing_igv hide variants that are missing IGV snapshot (default: FALSE)
#' @param pubmed_keywords optional variable of keywords to add to gene symbol for automatically generated pubmed search links
# #' @examples
# #' ***TODO***

create_candidate_table_html <- function(candidates, output_dir, sampleID, GTEx_median_rpkm=NULL, hide_missing_igv=FALSE, pubmed_keywords=NULL)
{
    candidates <- as.data.frame(candidates)
    if (hide_missing_igv) {
        candidates <- candidates[file.exists(candidates$igv_filename), ]
    }

    if (!dir.exists(paste0(output_dir, "html_files/candidate_variants/images/"))) {
        dir.create(paste0(output_dir, "html_files/candidate_variants/images/"), recursive=TRUE)
    }

    # Add additional columns or reformat for display output
    candidates[, "chr"] <- gsub("chr", "", candidates[, "chromosome"])
    candidates[, "base"] <- candidates[, "start"]
    candidates[, "ref"] <- candidates[, "reference"]
    candidates[, "alt"] <- candidates[, "alternate"]
    candidates$change <- gsub(".", "", gsub("non", "N", gsub("frameshift", "FS", gsub("synonymous", "S", candidates$change))), fixed=TRUE)
    candidates$change <- gsub("ertion", "", gsub("etion", "", candidates$change))
    candidates$SIFT[candidates$SIFT %in% c(".", "n/a", "not scored")] <- ""
    candidates$SIFT <- gsub(" *warning! low confidence.", "*", gsub(" due to stop", "", candidates$SIFT), fixed=TRUE)
    candidates$SIFT <- gsub("damaging", "D", gsub("tolerated", "T", candidates$SIFT))
    candidates$PP2 <- gsub(".", "", gsub("possibly damaging", "P", gsub("probably damaging", "D", gsub("benign", "B", candidates$Polyphen2))), fixed=TRUE)
    candidates$RVIS <- round(candidates$RVIS, 1)
    # *** Must be a cleaner way to do the below... ***
    candidates[, "MAF ExAC"] <- gsub("e-0", "e-", gsub("0.0e+00", "0", format(candidates[, "MAF ExAC"], digits=2), fixed=TRUE))
    candidates[, "MAF 1000G"] <- gsub("e-0", "e-", gsub("0.0e+00", "0", format(candidates[, "MAF 1000G"], digits=2), fixed=TRUE))
    
    output_cols <- c("inheritance model", "gene", "chr", "ref", "alt", "region", "change","ExAC count", 
                     "MAF ExAC", "MAF 1000G", "SIFT", "PP2", "RVIS")
    for (id in sampleID) {
        output_cols <- c(output_cols, paste(id, "genotype"))
    }
    table_candidates <- as.data.frame(candidates[, output_cols])

    # Abbreviate long sequences with "..." and hovertext
    table_candidates$ref <- sapply(table_candidates$ref, function(x){ifelse(nchar(x) > 5, paste0('<span title="', x, '">', substr(x, 1, 3), "... </span>"), x)})
    table_candidates$alt <- sapply(table_candidates$alt, function(x){ifelse(nchar(x) > 5, paste0('<span title="', x, '">', substr(x, 1, 3), "... </span>"), x)})

    # Convert columns to factors for nice filtering output
    table_candidates$"inheritance model" <- factor(table_candidates$"inheritance model")
    table_candidates$gene <- factor(table_candidates$gene)
    table_candidates$chr <- factor(table_candidates$chr)
    table_candidates$ref <- factor(table_candidates$ref)
    table_candidates$alt <- factor(table_candidates$alt)
    table_candidates$region <- factor(table_candidates$region)
    table_candidates$change <- factor(table_candidates$change)
    table_candidates[, "ExAC count"] <- as.numeric(table_candidates[, "ExAC count"])
    table_candidates[, "MAF ExAC"] <- as.numeric(table_candidates[, "MAF ExAC"])
    table_candidates[, "MAF 1000G"] <- as.numeric(table_candidates[, "MAF 1000G"])
    table_candidates$SIFT <- factor(table_candidates$SIFT)
    table_candidates$PP2 <- factor(table_candidates$PP2)
    table_candidates[, "RVIS"] <- as.numeric(table_candidates[, "RVIS"])
    for (id in sampleID) {
        id_colname <- paste(id, "genotype")
        table_candidates[, id_colname] <- factor(table_candidates[, id_colname])
    }

    # # *** ADD FIRST COLUMN "view" with links to individual candidates and FINAL COLUMN WITH EXTERNAL LINKS ***
    # view_col <- paste0("<a href='candidate_variants/candidate", 1:nrow(table_candidates), ".html'>view</a>")
    # hide_col <- rep("<input type='checkbox' id='checkbox'>", nrow(table_candidates))
    # table_candidates <- cbind(table_candidates[, 1], view_col, hide_col, table_candidates[, -1])
    # colnames(table_candidates)[2:3] <- c("view", "hide")

    rownames(table_candidates) <- NULL
    colnames(table_candidates) <- sapply(colnames(table_candidates), function(x) {paste(strwrap(x, width=5),  collapse="\n")})
    colnames(table_candidates) <- gsub("genotype", "GT", colnames(table_candidates))

    # IGV plot
    IGV_col <- ifelse(file.exists(candidates$igv_filename), paste0("<a href='candidate_variants/igv_output/", basename(candidates$igv_filename), "'>IGV</a>"), "")
    table_candidates <- cbind(table_candidates, IGV_col)
    colnames(table_candidates)[ncol(table_candidates)] <- "IGV\nsnapshot"

    # GTEx plot
    for (gg in unique(table_candidates$gene)) {
        gg_filename <- paste0(output_dir, "html_files/candidate_variants/images/GTEx_", gg, ".png")

        if (!file.exists(gg_filename)) {
            gg_p <- plot_gtex_expression(gg, GTEx_median_rpkm=GTEx_median_rpkm)
            if (!is.null(gg_p)) {
                ggplot2::ggsave(filename=gg_filename, plot=gg_p)
            }
        }
    }
    GTEx_col_filenames <- paste0(output_dir, "html_files/candidate_variants/images/GTEx_", table_candidates$gene, ".png")
    GTEx_col <- ifelse(file.exists(GTEx_col_filenames), paste0("<a href='candidate_variants/images/GTEx_", table_candidates$gene, ".png'>GTEx plot</a>"), "")
    table_candidates <- cbind(table_candidates, GTEx_col)
    colnames(table_candidates)[ncol(table_candidates)] <- "GTEx\nbarplot"

    # dbSNP link
    dbSNP_link <- paste0("<a href='https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?searchType=adhoc_search&type=rs&rs=", candidates$dbSNP, "'>", candidates$dbSNP, "</a>")
    dbSNP_col <- ifelse(candidates$dbSNP == ".", "", dbSNP_link)
    table_candidates <- cbind(table_candidates, dbSNP_col)
    colnames(table_candidates)[ncol(table_candidates)] <- "dbSNP"

    # ExAC link
    ExAC_SNV_link <- paste0("<a href='http://exac.broadinstitute.org/variant/", gsub("chr", "", candidates$chr), "-", candidates$start, "-", candidates$ref, "-", candidates$alt, "'>ExAC</a>")
    ExAC_region_link <- paste0("<a href='http://exac.broadinstitute.org/region/", gsub("chr", "", candidates$chr), "-", candidates$start, "-", candidates$start, "'>ExAC</a>")
    ExAC_link <- ifelse(endsWith(as.character(table_candidates$change), "SNV"), ExAC_SNV_link, ExAC_region_link)
    
    # gnomAD link (copied from ExAC link)
    gnomAD_link <- gsub("ExAC", "gnomAD", gsub("exac", "gnomad", ExAC_link))

    # OMIM link
    OMIM_mimid <- rep("", nrow(table_candidates))
    for (ii in 1:nrow(table_candidates)) {
        gg <- table_candidates$gene[ii]
        if (gg %in% mim2gene[, 4]) {
            gg_OMIM <- mim2gene[mim2gene[, 4] == gg, 1]
            if (length(gg_OMIM) == 1) {
                OMIM_mimid[ii] <- gg_OMIM
            }

        }
    }
    OMIM_link <- ifelse(OMIM_mimid == "", "", paste0("<a href='http://omim.org/entry/", OMIM_mimid, "'>OMIM</a>"))

    # GWAS Central
    GWAScentral_link <- paste0("<a href='http://www.gwascentral.org/generegion?q=", candidates$gene, "'>GWAScentral</a>")

    # Genecards
    GeneCards_link <- paste0("<a href='http://www.genecards.org/cgi-bin/carddisp.pl?gene=", candidates$gene, "'>GeneCards</a>")

    # Malacards
    Malacards_link <- paste0("<a href='http://www.malacards.org/search/results/", candidates$gene, "'>Malacards</a>")

    # # NCBI gene  *** TODO: Fix
    # symDB <- as.data.frame(org.Hs.eg.db::org.Hs.egSYMBOL)
    # Entrez_geneid <- symDB$gene_id[match(candidates$gene, symDB$symbol)]
    # NCBIgene_link <- ifelse(is.na(Entrez_geneid), "", paste0("<a href='https://www.ncbi.nlm.nih.gov/gene/", Entrez_geneid, "'>NCBI Gene</a>"))

    # Pubmed *** TODO: make this a link to gene entry instead of just a pubmed search for gene symbol ***
    pubmed_search_terms <- ifelse(is.null(pubmed_keywords), "", paste0("+", paste0(pubmed_keywords, collapse="+")))
    pubmed_link <- paste0("<a href='https://www.ncbi.nlm.nih.gov/pubmed?EntrezSystem2.PEntrez.Pubmed.SearchBar.Db=pubmed&term=", candidates$gene, pubmed_search_terms, "'>pubmed search: ", candidates$gene, pubmed_search_terms, "</a>")

    external_links <- paste(ExAC_link, gnomAD_link, OMIM_link, GWAScentral_link, GeneCards_link, Malacards_link, pubmed_link, sep=" <br> ")
    # *** FIND BETTER WAY TO DO BELOW... JUST UGLY HACK TO BRUTE FORCE IT FOR NOW ***
    external_links <- gsub(" <br>  <br> ", " <br> ", external_links, fixed=TRUE)
    external_links <- gsub(" <br>  <br> ", " <br> ", external_links, fixed=TRUE)
    external_links <- gsub(" <br>  <br> ", " <br> ", external_links, fixed=TRUE)
    external_links <- gsub(" <br>  <br> ", " <br> ", external_links, fixed=TRUE)
    external_links <- gsub(" <br>  <br> ", " <br> ", external_links, fixed=TRUE)
    external_links <- gsub(" <br>  <br> ", " <br> ", external_links, fixed=TRUE)
    external_links <- gsub(" <br>  <br> ", " <br> ", external_links, fixed=TRUE)
    external_links <- gsub(" <br>  <br> ", " <br> ", external_links, fixed=TRUE)
    table_candidates <- cbind(table_candidates, external_links)
    colnames(table_candidates)[ncol(table_candidates)] <- "External Links"

    candidate_types <- unique(as.character(table_candidates[, "inheritance\nmodel"]))
    for (ct in candidate_types) {
        filename <- paste0(output_dir, "html_files/candidate_variants_", gsub(" ", "_", ct), ".html")
        ct_table_candidates <- table_candidates[as.character(table_candidates[, "inheritance\nmodel"]) %in% ct, ]

        dt <- DT::datatable(ct_table_candidates, filter='bottom', rownames=FALSE, escape=FALSE, 
                                        extensions=c('Buttons'), #, 'Responsive'),
                                        options=list(dom='Blfrtip', buttons=c(I('colvis'), 'csv')))
        # *** TODO: FIX and EXPAND
        # dt_format <- dt %>% 
        #             formatStyle("SIFT", color=styleEqual(c("T", "D", "D*", "."), c("green", "red", "red", "black"))) %>%
        #             formatStyle("PP2", color=styleEqual(c("B", "P", "D", "."), c("green", "orange", "red", "black")))

        DT::saveWidget(dt, filename, selfcontained=FALSE)
    }
}

