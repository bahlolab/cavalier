
#' @importFrom readr read_delim cols
#' @importFrom dplyr "%>%" mutate rename
get_gtex_expression <- function()
{
    gtex_gene_median_tpm <- cavalier_cache$gtex_gene_median_tpm
    
    if (is.null(gtex_gene_median_tpm)) {
        
        gtex_gene_median_tpm_uri <- get_cavalier_opt('gtex_gene_median_tpm_uri')
        cache_dir <- get_cavalier_opt('cache_dir')
        fn <- basename(gtex_gene_median_tpm_uri) %>% 
            str_remove('.txt$') %>% 
            str_c('.rds') %>% 
            file.path(cache_dir, .)
        
        gtex_gene_median_tpm <- 
            (function() read_delim(gtex_gene_median_tpm_uri,
                                   delim = "\t",
                                   skip = 2,
                                   col_types = cols())) %>% 
            cache(fn) %>% 
            rename(ensembl_gene_id = Name,
                   symbol = Description) %>% 
            mutate(ensembl_gene_id = str_remove(ensembl_gene_id, '\\.[0-9]+$'),
                   symbol = {
                       s1 <- hgnc_ensembl2sym(ensembl_gene_id)
                       at <- is.na(s1)
                       replace(s1, at, hgnc_sym2sym(symbol[at], remove_unknown = TRUE))})

        cavalier_cache$gtex_gene_median_tpm <- gtex_gene_median_tpm
    }
    
    return(gtex_gene_median_tpm)
}

#' @export
get_gtex_tissues <- function()
{
    get_gtex_expression() %>% 
        select(-(1:2)) %>% 
        colnames()
}

#' Plot GTEx tissue median RPKM expression for given gene symbol
#' 
#' @param gene gene symbol
#' @param gtex_gene_median_tpm location of GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct.gz file downloaded from GTEx Portal (see https://gtexportal.org/home/datasets)
#' @param small_font TRUE or FALSE to use smaller font size (for PDF output instead of HTML)
#' @param tissues optionally specify list of tissues to plot
#' @return ggplot2 plot of median RPKM expression for gene
# #' @examples
# #' ***TODO***
plot_gtex_expression <- function(gene, gtex_gene_median_tpm,
                                 small_font = FALSE,
                                 tissues = NULL)
{

    # Select tissue types to use
    if (is.null(tissues)) {
        # Specify list of default tissues to use
        tissues <- c(
                        "Brain - Amygdala",
                        "Brain - Anterior cingulate cortex (BA24)",
                        "Brain - Caudate (basal ganglia)",
                        "Brain - Cerebellar Hemisphere",
                        "Brain - Cerebellum",
                        "Brain - Cortex",
                        "Brain - Frontal Cortex (BA9)",
                        "Brain - Hippocampus",
                        "Brain - Hypothalamus",
                        "Brain - Nucleus accumbens (basal ganglia)",
                        "Brain - Putamen (basal ganglia)",
                        "Brain - Spinal cord (cervical c-1)",
                        "Brain - Substantia nigra",
                        "Adipose - Subcutaneous",
                        "Artery - Tibial",
                        "Breast - Mammary Tissue",
                        "Cells - Transformed fibroblasts",
                        "Esophagus - Mucosa",
                        "Lung",
                        "Muscle - Skeletal",
                        "Nerve - Tibial",
                        "Skin - Not Sun Exposed (Suprapubic)",
                        "Skin - Sun Exposed (Lower leg)",
                        "Thyroid", 
                        "Whole Blood"
                    )
        
        # # Alternative: use all tissue types
        # tissues <- colnames(GTEx_full)[10:ncol(GTEx_full)]
        # # Sort columns to put brain tissues first
        # tissues <- c(tissues[substr(tissues, 1, 5) == "Brain"], tissues[substr(tissues, 1, 5) != "Brain"])
    }

    # Match gene
    expression <- gtex_gene_median_tpm_table[gtex_gene_median_tpm_table$Description == gene, tissues]
    if (nrow(expression) == 0) {
        print(paste("Warning: zero rows found for gene:", gene))
        print("No plot returned")
        return(NULL)
    } else if (nrow(expression) > 1) {
        print(paste("Warning:", nrow(expression), "rows found for gene:", gene))
        print("Plotting mean expression across all rows")
        expression <- colMeans(expression)
    }
    
    # Plot log2rpkm expression of tissue types for gene
    expression_df <- data.frame(tissue=tissues, expression=log2(as.numeric(as.vector(unlist(expression)))+1), stringsAsFactors=FALSE)
    expression_df$tissue <- factor(expression_df$tissue, levels=tissues)
    expression_df$tissue_name <- gsub(" - ", ": ", as.character(expression_df$tissue))
    expression_df$tissue_name <- factor(expression_df$tissue_name, levels=gsub(" - ", ": ", tissues))
    expression_df$tissue_class <- sapply(as.character(expression_df$tissue), function(x){strsplit(x, " - ")[[1]][1]})

    p <- 
        ggplot2::ggplot(data=expression_df, ggplot2::aes(x=tissue_name, y=expression, fill=tissue_class)) + 
        ggplot2::geom_bar(stat="identity") + 
        ggplot2::scale_fill_manual(values=c(RColorBrewer::brewer.pal(8, "Set1"), RColorBrewer::brewer.pal(8, "Dark2"))[c(1, 3, 2, 4:16)]) + 
        ggplot2::ggtitle(paste(gene, "GTEx")) + 
        ggplot2::ylab("log-2 RPKM") + ggplot2::xlab("") +
        ggplot2::theme_bw() + 
        ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, hjust=1, vjust=0.25)) + 
        ggplot2::guides(fill=FALSE) +
        ggplot2::coord_flip()
    if (small_font) {
        p <- p + ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, hjust=1, vjust=0.25, size=6), axis.text.y=ggplot2::element_text(size=6), axis.title.y=ggplot2::element_text(size=6), title=ggplot2::element_text(size=8))
    }
    return(p)
}


