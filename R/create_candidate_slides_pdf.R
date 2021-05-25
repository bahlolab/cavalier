#' Create PDF output slide for each candidate variant
#' 
#' @param candidates candidate variants data.frame
#' @param output_dir cavalier output directory
#' @param genemap2 location of genemap2.txt file downloaded from OMIM (see https://omim.org/downloads/)
#' @param GTEx_median_rpkm location of GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct.gz file downloaded from GTEx Portal (see https://gtexportal.org/home/datasets)
#' @param GTEx_tissues optionally specify list of tissues to plot GTEx expression data
#' @param hide_missing_igv hide variants that are missing IGV snapshot (default: FALSE)
#' @param layout slide layout choice: "individual" or "multiple" designed for a single or multiple individuals (default: "individual")
# #' @examples
# #' ***TODO***

create_candidate_slides_pdf <- function(candidates, output_dir, output_cols, genemap2=NULL, GTEx_median_rpkm=NULL, GTEx_tissues=NULL, hide_missing_igv=FALSE, layout="individual")
{
    # If option specified then remove candidates without IGV files
    if (hide_missing_igv) {
        candidates <- candidates[file.exists(candidates$igv_filename), ]
    }
    if (nrow(candidates) == 0) {
        return(NULL)
    }
    rownames(candidates) <- 1:nrow(candidates)
    
    if (all(c("reference", "alternate") %in% colnames(candidates))) {
        candidates$"ref / alt" <- paste(candidates$reference, "/", candidates$alternate)
    }
    if ("change" %in% colnames(candidates)) {
        candidates$change <- gsub(";", "\n", candidates$change, fixed=TRUE)
    }
    
    output_cols <- intersect(output_cols, colnames(candidates))
    
    output_dir <- endslash_dirname(output_dir)
    
    if (!dir.exists(paste0(output_dir, "pdf_files/"))) {
        dir.create(paste0(output_dir, "pdf_files/"), recursive=TRUE)
    }
    
    # Format strings for multiple lines
    max_line_length <- 20
    candidates$"ref / alt" <- sapply(candidates$"ref / alt", function(x){newline_every_n_chars(x, max_line_length)})

    table_base_size <- 14
    omim_table_base_size <- 12
    
    candidate_types <- unique(as.character(candidates[, "inheritance model"]))
    for (ct in candidate_types) {
        filename <- paste0(output_dir, "pdf_files/candidate_variants_report_", gsub(" ", "_", ct), ".pdf")
        ct_candidates <- candidates[as.character(candidates$`inheritance model`) %in% ct, ]

        pdf(filename, 10, 8)
        for (ii in 1:nrow(ct_candidates)) {
            ii_gene <- ct_candidates$gene[ii]

            # Title
            ii_title <- paste0("Candidate ", ii, ":  ", ii_gene)
            ii_title_grob <- grid::textGrob(ii_title, gp=grid::gpar(cex=2, col="navy"))
            
            # Information table
            ii_table <- data.frame(value=as.vector(unlist(ct_candidates[ii, output_cols])))
            rownames(ii_table) <- output_cols
            ii_table_grob <- gridExtra::tableGrob(ii_table, cols=NULL, theme=gridExtra::ttheme_minimal(base_size=table_base_size, 
                                                 core=list(fg_params=list(hjust=0, x=0.05))))
            # *** TODO: colour text for different predictions green / red etc
            
            # IGV
            if (file.exists(ct_candidates$igv_filename[ii])) {
                ii_igv_png <- png::readPNG(ct_candidates$igv_filename[ii], native=TRUE, info=TRUE)
                ii_igv_grob <- grid::rasterGrob(ii_igv_png)
            } else {
                ii_igv_grob <- grid::textGrob("")
            }

            # GTEx
            ii_GTEx <- plot_gtex_expression(ii_gene, GTEx_median_rpkm=GTEx_median_rpkm, small_font=TRUE, tissues=GTEx_tissues)
            if (is.null(ii_GTEx)) {
                ii_GTEx <- grid::textGrob("")
            }

            # *** OMIM ***
            ii_omim_table <- omim_table(ii_gene, genemap2=genemap2)
            if (is.null(ii_omim_table)) {
                ii_omim_grob <- grid::textGrob("")
            }
            else {
                ii_omim_grob <- gridExtra::tableGrob(ii_omim_table, rows=NULL,
                                    theme=gridExtra::ttheme_default(base_size=omim_table_base_size, 
                                                         core=list(fg_params=list(hjust=0, x=0.05))))
            }
            
            # Combine images

            # # *** DEFINE SEVERAL LAYOUTS ***
            # if (length(project_info$sampleID) == 1) {
            #     grid_layout <- rbind(c(1, 1, 1),
            #                          c(2, 3, 3),
            #                          c(4, 4, 5))
            # } 
            # else {
            #     grid_layout <- rbind(c(1, 1, 1),
            #                          c(2, 5, 3),
            #                          c(4, 4, 3))
            # }
            # grid_widths <- c(4, 2, 3)
            # grid_heights <- c(1, 5, 3)

            # *** DEFINE SEVERAL LAYOUTS ***
            if (layout == "individual") {
                grid_layout <- rbind(c(6,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  7),
                                     c(6,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  7),
                                     c(6,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  7),
                                     c(6,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  7),
                                     c(6,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9, 12,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  7),
                                     c(6,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4, 12,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  7),
                                     c(6,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4, 12,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  7),
                                     c(6,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4, 12,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  7),
                                     c(6,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4, 12,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  7),
                                     c(6,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4, 12,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  7),
                                     c(6,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4, 12,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  7),
                                     c(6,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4, 12,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  7),
                                     c(6,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4, 12,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  7),
                                     c(6,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  7))
            } 
            else if (layout == "multiple") {
                grid_layout <- rbind(c(6,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  7),
                                     c(6,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  7),
                                     c(6,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  7),
                                     c(6,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10,  5,  5,  5,  5,  5,  5,  5,  5, 13,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10,  5,  5,  5,  5,  5,  5,  5,  5, 13,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10,  5,  5,  5,  5,  5,  5,  5,  5, 13,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10,  5,  5,  5,  5,  5,  5,  5,  5, 13,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10,  5,  5,  5,  5,  5,  5,  5,  5, 13,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10,  5,  5,  5,  5,  5,  5,  5,  5, 13,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10,  5,  5,  5,  5,  5,  5,  5,  5, 13,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10,  5,  5,  5,  5,  5,  5,  5,  5, 13,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10,  5,  5,  5,  5,  5,  5,  5,  5, 13,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10,  5,  5,  5,  5,  5,  5,  5,  5, 13,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10,  5,  5,  5,  5,  5,  5,  5,  5, 13,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10,  5,  5,  5,  5,  5,  5,  5,  5, 13,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10,  5,  5,  5,  5,  5,  5,  5,  5, 13,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10,  5,  5,  5,  5,  5,  5,  5,  5, 13,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10,  5,  5,  5,  5,  5,  5,  5,  5, 13,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9, 12,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4, 12,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4, 12,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4, 12,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4, 12,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4, 12,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4, 12,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4, 12,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4, 12,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  7),
                                     c(6,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  7))
            } else {
                stop("Error: unrecognised layout option")
            }
            grid_widths <- rep(0.25, 40)
            grid_heights <- rep(0.25, 32)

            ii_grob_list <- list(ii_title_grob, ii_table_grob, ii_igv_grob, ii_omim_grob, ii_GTEx, grid::textGrob(""), grid::textGrob(""), grid::textGrob(""), grid::textGrob(""), grid::textGrob(""), grid::textGrob(""), grid::textGrob(""), grid::textGrob(""))
            gridExtra::grid.arrange(grobs=ii_grob_list, layout_matrix=grid_layout, widths=grid_widths, heights=grid_heights, right="")

            grid::textGrob("")
            grid::textGrob("")
        }
        dev.off()
    }
}

