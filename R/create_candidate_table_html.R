#' Create HTML table of candidate variants
#' 
#' @param candidates candidate variants data.frame
#' @param output_dir cavalier output directory
#' @param sampleID list of IDs and names for samples of interest
#' @param GTEx_median_rpkm location of GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct.gz file downloaded from GTEx Portal (see https://gtexportal.org/home/datasets)
#' @param GTEx_tissues optionally specify list of tissues to plot GTEx expression data
#' @param hide_missing_igv hide variants that are missing IGV snapshot (default: FALSE)
#' @param pubmed_keywords optional variable of keywords to add to gene symbol for automatically generated pubmed search links
# #' @examples
# #' ***TODO***

create_candidate_table_html <- function(candidates, output_dir, sampleID, output_cols, GTEx_median_rpkm=NULL, GTEx_tissues=NULL, hide_missing_igv=FALSE, pubmed_keywords=NULL)
{
    # If option specified then remove candidates without IGV files
    if (hide_missing_igv) {
        candidates <- candidates[file.exists(candidates$igv_filename), ]
    }
    if (nrow(candidates) == 0) {
        return(NULL)
    }
    
    output_dir <- endslash_dirname(output_dir)

    if (!dir.exists(paste0(output_dir, "html_files/candidate_variants/images/"))) {
        dir.create(paste0(output_dir, "html_files/candidate_variants/images/"), recursive=TRUE)
    }
    
    output_cols <- unique(c(intersect(output_cols, colnames(candidates)), "inheritance model"))
    
    for (id in sampleID) {
        output_cols <- c(output_cols, paste(id, "genotype"))
    }
    table_candidates <- as.data.frame(candidates[, output_cols])
    
    # *** GENERALISE
    # # Abbreviate long sequences with "..." and hovertext
    # table_candidates$ref <- sapply(table_candidates$ref, function(x){ifelse(nchar(x) > 5, paste0('<span title="', x, '">', substr(x, 1, 3), "... </span>"), x)})
    # table_candidates$alt <- sapply(table_candidates$alt, function(x){ifelse(nchar(x) > 5, paste0('<span title="', x, '">', substr(x, 1, 3), "... </span>"), x)})

    # Convert columns to factors for nice filtering output
    table_candidates$"inheritance model" <- factor(table_candidates$"inheritance model")
    for (id in sampleID) {
        id_colname <- paste(id, "genotype")
        table_candidates[, id_colname] <- factor(table_candidates[, id_colname])
    }
    
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
            gg_p <- plot_gtex_expression(gg, GTEx_median_rpkm=GTEx_median_rpkm, tissues=GTEx_tissues)
            if (!is.null(gg_p)) {
                ggplot2::ggsave(filename=gg_filename, plot=gg_p)
            }
        }
    }
    GTEx_col_filenames <- paste0(output_dir, "html_files/candidate_variants/images/GTEx_", table_candidates$gene, ".png")
    GTEx_col <- ifelse(file.exists(GTEx_col_filenames), paste0("<a href='candidate_variants/images/GTEx_", table_candidates$gene, ".png'>GTEx plot</a>"), "")
    table_candidates <- cbind(table_candidates, GTEx_col)
    colnames(table_candidates)[ncol(table_candidates)] <- "GTEx\nbarplot"

    # gnomAD link
    gnomAD_SNV_link <- paste0("<a href='http://gnomAD.broadinstitute.org/variant/", gsub("chr", "", candidates$chr), "-", candidates$start, "-", candidates$ref, "-", candidates$alt, "'>gnomAD</a>")
    gnomAD_region_link <- paste0("<a href='http://gnomAD.broadinstitute.org/region/", gsub("chr", "", candidates$chr), "-", candidates$start, "-", candidates$start, "'>gnomAD</a>")
    gnomAD_link <- ifelse(nchar(candidates$ref) == 1 & nchar(candidates$alt) == 1, gnomAD_SNV_link, gnomAD_region_link)
    
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

    # Genecards
    GeneCards_link <- paste0("<a href='http://www.genecards.org/cgi-bin/carddisp.pl?gene=", candidates$gene, "'>GeneCards</a>")
    
    # Mouse Genome Informatics
    MGI_link <- paste0("<a href='http://www.informatics.jax.org/searchtool/Search.do?query=", candidates$gene, "&submit=Quick%0D%0ASearch'>MGI</a>")

    # Pubmed *** TODO: make this a link to gene entry instead of just a pubmed search for gene symbol ***
    pubmed_search_terms <- ifelse(is.null(pubmed_keywords), "", paste0("+", paste0(pubmed_keywords, collapse="+")))
    pubmed_link <- paste0("<a href='https://www.ncbi.nlm.nih.gov/pubmed?EntrezSystem2.PEntrez.Pubmed.SearchBar.Db=pubmed&term=", candidates$gene, pubmed_search_terms, "'>pubmed search: ", candidates$gene, pubmed_search_terms, "</a>")

    external_links <- paste(gnomAD_link, OMIM_link, GeneCards_link, MGI_link, pubmed_link, sep=" <br> ")
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
        if (substring(filename, 1, 1) != "/") {
            filename <- paste0(getwd(), "/", filename)
        }
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

