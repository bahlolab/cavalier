#' Create PDF slide(s) containing table of candidate variants
#' 
#' @param candidates candidate variants data.frame
#' @param output_dir cavalier output directory
#' @param hide_missing_igv hide variants that are missing IGV snapshot (default: FALSE)
# #' @examples
# #' ***TODO***

create_candidate_table_pdf <- function(candidates, output_dir, output_cols, hide_missing_igv=FALSE)
{
    # If option specified then remove candidates without IGV files
    if (hide_missing_igv) {
        candidates <- candidates[file.exists(candidates$igv_filename), ]
    }
    if (nrow(candidates) == 0) {
        return(NULL)
    }
    
    output_dir <- endslash_dirname(output_dir)
    
    if (!dir.exists(paste0(output_dir, "pdf_files/"))) {
        dir.create(paste0(output_dir, "pdf_files/"), recursive=TRUE)
    }
    
    candidates$chr <- candidates$chromosome
    candidates$ref <- sapply(candidates$reference, function(x){ifelse(nchar(x) > 7, paste0(substr(x, 1, 5), "..."), x)})
    candidates$alt <- sapply(candidates$alternate, function(x){ifelse(nchar(x) > 7, paste0(substr(x, 1, 5), "..."), x)})
    
    output_cols <- intersect(output_cols, colnames(candidates))

    table_candidates <- candidates[, output_cols]
    rownames(table_candidates) <- NULL

    if (!dir.exists(paste0(output_dir, "pdf_files/"))) {
        dir.create(paste0(output_dir, "pdf_files/"), recursive=TRUE)
    }

    candidate_types <- unique(candidates[, "inheritance model"])
    for (ct in candidate_types) {
        title <- paste(ct, "candidate variants")
        filename <- paste0(output_dir, "pdf_files/candidate_variants_table_", gsub(" ", "_", ct), ".pdf")
        ct_table_candidates <- table_candidates[candidates[, "inheritance model"] %in% ct, ]
        rownames(ct_table_candidates) <- NULL

        pdf(filename, width=10, height=8)
        if (nrow(ct_table_candidates) <= 20) {
            title_grob <- grid::textGrob(title, gp=grid::gpar(cex=2, col="navy"))
            table_grob <- gridExtra::tableGrob(ct_table_candidates, theme=gridExtra::ttheme_default(base_size=12))

            gridExtra::grid.arrange(grobs=list(title_grob, table_grob), layout_matrix=rbind(1, 2), heights=c(1, 8))
            }
        else {
            rows_per_page <- 20
            pages <- ceiling(nrow(ct_table_candidates) / rows_per_page)
            
            for (pp in 1:pages) {
                pp_title_grob <- grid::textGrob(paste0(title, " (", pp, " of ", pages, ")"), gp=grid::gpar(cex=2, col="navy"))
                pp_start_row <- 1 + rows_per_page * (pp - 1)
                pp_table_rows <- pp_start_row:min((pp_start_row + (rows_per_page - 1)), nrow(ct_table_candidates))
                pp_table_grob <- gridExtra::tableGrob(ct_table_candidates[pp_table_rows, ], theme=gridExtra::ttheme_default(base_size=12))
                gridExtra::grid.arrange(grobs=list(pp_title_grob, pp_table_grob), layout_matrix=rbind(1, 2), heights=c(1, 8))
            }
        }
        dev.off()
    }
}

