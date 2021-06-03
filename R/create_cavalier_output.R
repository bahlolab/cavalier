#' Create cavalier output
#' 
#' @param candidates candidate variants data.frame
#' @param output_dir cavalier output directory
#' @param sampleID list of IDs and names for samples of interest
#' @param output_cols list of columns to output
#' @param hide_missing_igv hide variants that are missing IGV snapshot (default: FALSE)
#' @param pubmed_keywords optional variable of keywords to add to gene symbol for automatically generated pubmed search links
#' @param layout slide layout choice: "individual" or "multiple" designed for a single or multiple individuals (default: "individual")
#' @param GTEx_median_rpkm location of GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct.gz file downloaded from GTEx Portal (see https://gtexportal.org/home/datasets)
#' @param GTEx_tissues optionally specify list of tissues to plot GTEx expression data
#' @param genemap2 location of genemap2.txt file downloaded from OMIM (see https://omim.org/downloads/)
# #' @examples
# #' ***TODO***

create_cavalier_output <- function(candidates, output_dir, sampleID, output_cols, hide_missing_igv=FALSE, pubmed_keywords="", layout="individual", genemap2=NULL, GTEx_median_rpkm=NULL, GTEx_tissues=NULL,
                                   add_data_col = NULL, title_col = NULL)
{
    output_dir <- endslash_dirname(output_dir)
    
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive=TRUE)
    }
    
    # Create output directory for HTML and PDF files if they do not exist
    if (!dir.exists(paste0(output_dir, "html_files/candidate_variants/"))) {
        dir.create(paste0(output_dir, "html_files/candidate_variants/"), recursive=TRUE)
    }
    if (!dir.exists(paste0(output_dir, "pdf_files/"))) {
        dir.create(paste0(output_dir, "pdf_files/"))
    }
    
    # Copy peddy.html and IGV plots to HTML output dir (if they exist)
    if (file.exists(paste0(output_dir, "data/peddy/peddy.html"))) {
        system(paste0("cp ", output_dir, "data/peddy/peddy.html ", output_dir, "/html_files/peddy.html"))
    }
    if (dir.exists(paste0(output_dir, "data/igv_output/"))) {
        # Remove existing IGV output dir (if it exists)
        if (dir.exists(paste0(output_dir, "html_files/candidate_variants/igv_output/"))) {
            system(paste0("rm -r ", output_dir, "html_files/candidate_variants/igv_output"))
        }
        system(paste0("cp -r ", output_dir, "data/igv_output/ ", output_dir, "html_files/candidate_variants/igv_output/"))
    }
    
    create_candidate_table_html(candidates, output_dir, sampleID, output_cols, hide_missing_igv=hide_missing_igv, 
                                pubmed_keywords=pubmed_keywords, GTEx_median_rpkm=GTEx_median_rpkm, GTEx_tissues=GTEx_tissues)

    create_candidate_table_pdf(candidates, output_dir, output_cols, hide_missing_igv=hide_missing_igv)

    create_candidate_slides_pdf(candidates, output_dir, output_cols, hide_missing_igv=hide_missing_igv, 
                                GTEx_median_rpkm=GTEx_median_rpkm, GTEx_tissues=GTEx_tissues, genemap2=genemap2, layout=layout,
                                add_data_col = add_data_col, title_col = title_col)
    
    # Write table of candidate variants
    if (hide_missing_igv) {
        candidates <- candidates[file.exists(candidates$igv_filename), ]
    }
    if (nrow(candidates) > 0) {
        if(!is.null(add_data_col)) {
            candidates[[add_data_col]] <- NULL
        }
        write.table(candidates, file=paste0(output_dir, "/candidate_variants.txt"), quote=FALSE, sep="\t", row.names=FALSE)
    }
}

