
#' @importFrom readr read_delim cols
#' @importFrom dplyr "%>%" mutate rename
get_gtex_expression <- function()
{
    gtex_gene_median_tpm <- cavalier_cache$gtex_gene_median_tpm
    
    if (is.null(gtex_gene_median_tpm)) {
        
        gtex_gene_median_tpm_uri <- get_cavalier_opt('gtex_gene_median_tpm_uri')
        cache_dir <- get_cavalier_opt('cache_dir')
        
        gtex_gene_median_tpm <- 
            (function() read_delim(gtex_gene_median_tpm_uri,
                                   delim = "\t",
                                   skip = 2,
                                   col_types = cols())) %>% 
            cache(basename(gtex_gene_median_tpm_uri)) %>% 
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

#' @importFrom cowplot ggdraw draw_text
#' @importFrom ggplot2 ggplot aes geom_bar scale_fill_manual ggtitle ylab xlab theme_bw theme guides coord_flip
#' @export
plot_gtex_expression <- function(gene, ensembl_id = NULL)
{
    tissues <- get_cavalier_opt('gtex_tissues')
    
    assert_that(is_scalar_character(gene),
                is.null(ensembl_id) || is_scalar_character(ensembl_id),
                is_character(tissues),
                all(tissues %in% get_gtex_tissues()))
    
    gtex_gene_median_tpm <- get_gtex_expression()
    
    # if gene not found return plot stating as such
    if ((is.null(ensembl_id) & !gene %in% gtex_gene_median_tpm$symbol) | 
        (!is.null(ensembl_id) & !ensembl_id %in% gtex_gene_median_tpm$ensembl_gene_id)) {
        return(ggdraw() + draw_text(str_c('"', gene, '"\n not found in GTEx')))
    }
    
    gene_exp <-
        `if`(!is.null(ensembl_id),
             gtex_gene_median_tpm %>% 
                 filter(ensembl_gene_id == ensembl_id) %>% 
                 select(all_of(tissues)),
             gtex_gene_median_tpm %>% 
                 filter(symbol == gene) %>% 
                 select(all_of(tissues))) %>% 
        mutate(id = seq_len(n())) %>% 
        pivot_longer(-id,
                     names_to = 'tissue',
                     values_to = 'expression') %>% 
        mutate(tissue = str_replace(tissue, ' -', ':') %>% str_remove('\\s?\\(.+\\)'),
               tissue_class = coalesce(str_extract(tissue, '^.+(?=:)'),
                                       tissue))
    
    if (n_distinct(gene_exp$id) == 1) {
        # plot single gene as barplot
        p <-
            gene_exp %>% 
            ggplot(aes(tissue, expression, fill=tissue_class)) + 
            geom_col() +
            labs(title = str_c("GTEx: ", gene),
                 y = "log-2 TPM") +
            theme_bw() + 
            theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.25),
                  axis.title.y = element_blank()) + 
            guides(fill=FALSE) +
            coord_flip()
    } else {
        # plot multiple genes as boxplot
        p <-
            gene_exp %>% 
            ggplot(aes(tissue, expression, col=tissue_class)) + 
            geom_boxplot() +
            labs(title = str_c("GTEx: ", gene),
                 y = "log-2 TPM") +
            theme_bw() + 
            theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.25),
                  axis.title.y = element_blank()) + 
            guides(fill=FALSE) +
            coord_flip()
    }
    
    return(p)
}


