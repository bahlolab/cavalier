#' Create cavalier output
#' 
#' @param candidates candidate variants data.frame
#' @param output_dir cavalier output directory
#' @param sampleID list of IDs and names for samples of interest
#' @param hide_missing_igv hide variants that are missing IGV snapshot (default: FALSE)
#' @param pubmed_keywords optional variable of keywords to add to gene symbol for automatically generated pubmed search links
#' @param layout slide layout choice: "individual" or "multiple" designed for a single or multiple individuals (default: "individual")
#' @param GTEx_median_rpkm location of GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct.gz file downloaded from GTEx Portal (see https://gtexportal.org/home/datasets)
#' @param genemap2 location of genemap2.txt file downloaded from OMIM (see https://omim.org/downloads/)
# #' @examples
# #' ***TODO***

create_cavalier_output <- function(candidates, output_dir, sampleID, hide_missing_igv=FALSE, pubmed_keywords="", layout="individual", genemap2=NULL, GTEx_median_rpkm=NULL)
{
    starting_wd <- getwd()
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive=TRUE)
    }
    setwd(output_dir)
    
    # Create output directory for HTML and PDF files if they do not exist
    if (!dir.exists("html_files/")) {
        dir.create("html_files/candidate_variants/", recursive=TRUE)
    }
    if (!dir.exists("pdf_files/")) {
        dir.create("pdf_files/")
    }
    
    # Copy peddy.html and IGV plots to HTML output dir (if they exist)
    if (file.exists("data/peddy/peddy.html")) {
        system("cp data/peddy/peddy.html html_files/peddy.html")
    }
    if (dir.exists("data/igv_output/")) {
        # Remove existing IGV output dir (if it exists)
        if (dir.exists("html_files/candidate_variants/igv_output/")) {
            system("rm -r html_files/candidate_variants/igv_output")
        }
        system("cp -r data/igv_output/ html_files/candidate_variants/igv_output/")
    }
    
    create_candidate_table_html(candidates, output_dir, sampleID, hide_missing_igv=hide_missing_igv, 
                                pubmed_keywords=pubmed_keywords, GTEx_median_rpkm=GTEx_median_rpkm)

    create_candidate_table_pdf(candidates, output_dir, hide_missing_igv=hide_missing_igv)

    create_candidate_slides_pdf(candidates, output_dir, hide_missing_igv=hide_missing_igv, 
                                GTEx_median_rpkm=GTEx_median_rpkm, genemap2=genemap2, layout=layout)
    
    # Write table of candidate variants
    if (hide_missing_igv) {
        write.table(candidates[file.exists(candidates$igv_filename), ], file="candvars.txt", quote=FALSE, sep="\t", row.names=FALSE)
    } else {
        write.table(candidates, file="candvars.txt", quote=FALSE, sep="\t", row.names=FALSE)
    }
    
    setwd(starting_wd)
}

