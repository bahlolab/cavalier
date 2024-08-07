
#' @importFrom readr read_delim cols
#' @importFrom dplyr "%>%" mutate rename
get_gtex_expression <- function()
{
  gtex_gene_median_tpm_url <- get_cavalier_opt("gtex_gene_median_tpm_url")
  
  stopifnot(is_scalar_character(gtex_gene_median_tpm_url))

  fun <- function() {
    parsed <- httr::parse_url(gtex_gene_median_tpm_url)
    
    if (parsed$scheme == 'file') {
      con <- parsed$path
    } else {
      con <-  
        retry('GET', parsed) %>% 
        content() %>% 
        rawConnection() %>% 
        gzcon()
    }
    read_delim(
      con,
      delim = "\t",
      skip = 2,
      col_types = cols()) %>% 
      rename(ensembl_gene_id = Name,
             symbol = Description) 
  }
  
  cache(
    fun = fun,
    name = 'GTEx_gene_median_tpm')

}

#' Replace symbol in GTEx expression with current HGNC version
get_gtex_expression_hgnc <- function()
{
  get_gtex_expression() %>% 
    mutate(ensembl_gene_id = str_remove(ensembl_gene_id, '\\.[0-9]+$'),
           symbol = coalesce(hgnc_ensembl2sym(ensembl_gene_id),
                             hgnc_sym2sym(symbol)))
  
}

#' @export
get_gtex_tissues <- function()
{
  get_gtex_expression_hgnc() %>% 
    select(-(1:2)) %>% 
    colnames()
}

#' Plot GTEx tissue median RPKM expression for given gene symbol
#' 

#' @importFrom cowplot ggdraw draw_text
#' @importFrom dplyr n_distinct
#' @importFrom ggplot2 ggplot aes geom_bar scale_fill_manual ggtitle ylab xlab theme_bw theme guides 
#' @importFrom ggplot2 coord_flip geom_col labs geom_text element_blank margin ylim facet_wrap element_text
#' @export
plot_gtex_expression <- function(gene, ensembl_id = NULL)
{
    tissues <- get_cavalier_opt('gtex_tissues')
    
    assert_that(is_scalar_character(gene),
                is.null(ensembl_id) || is_scalar_character(ensembl_id),
                is_character(tissues),
                all(tissues %in% get_gtex_tissues()))
    
    gtex_gene_median_tpm <- get_gtex_expression_hgnc()
    
    if (!is.null(ensembl_id) && 
        !ensembl_id %in% gtex_gene_median_tpm$ensembl_gene_id) {
        ensembl_id <- NULL
    }
    
    # if gene not found return plot stating as such
    if ((is.null(ensembl_id) && !gene %in% gtex_gene_median_tpm$symbol) | 
        (!is.null(ensembl_id) && !ensembl_id %in% gtex_gene_median_tpm$ensembl_gene_id)) {
      return(NULL)
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
            guides(fill='none') +
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
            guides(fill='none') +
            coord_flip()
    }
    
    return(p)
}

# Note: WIP
#' @importFrom cowplot ggdraw draw_text
#' @importFrom forcats as_factor
#' @importFrom ggplot2 ggplot aes geom_bar scale_fill_manual ggtitle ylab xlab theme_bw theme guides coord_flip
#' @importFrom cowplot plot_grid
#' @export
plot_gtex_compact <- function(gene, ensembl_id = NULL, top_n = 3)
{
    tissues <- get_cavalier_opt('gtex_tissues')
    
    assert_that(is_scalar_character(gene),
                is.null(ensembl_id) || is_scalar_character(ensembl_id),
                is_character(tissues),
                all(tissues %in% get_gtex_tissues()))
    
    gtex_gene_median_tpm <- get_gtex_expression_hgnc()
    
    # if gene not found return plot stating as such
    if ((is.null(ensembl_id) && !gene %in% gtex_gene_median_tpm$symbol) | 
        (!is.null(ensembl_id) && !ensembl_id %in% gtex_gene_median_tpm$ensembl_gene_id)) {
        return(ggdraw() + draw_text(str_c('"', gene, '"\n not found in GTEx')))
    }
    
    gene_exp <-
        `if`(!is.null(ensembl_id),
             gtex_gene_median_tpm %>% 
                 filter(ensembl_gene_id == ensembl_id) %>% 
                 select(all_of(get_gtex_tissues())),
             gtex_gene_median_tpm %>% 
                 filter(symbol == gene) %>% 
                 select(all_of(get_gtex_tissues()))) %>% 
        mutate(id = seq_len(n())) %>% 
        pivot_longer(-id,
                     names_to = 'tissue',
                     values_to = 'expression') %>% 
        mutate(tissue = str_replace(tissue, ' -', ':') %>% str_remove('\\s?\\(.+\\)'),
               tissue = coalesce(str_extract(tissue, '^.+(?=:)'),
                                       tissue)) %>% 
        group_by(tissue) %>% 
        summarise(expression = mean(expression),
                  .groups = 'drop') %>% 
        filter(tissue != 'Cells') %>% 
        arrange(expression) %>% 
        mutate(tissue = as_factor(tissue))
    
    ymx <- max(gene_exp$expression)
        
    
    # if (n_distinct(gene_exp$id) == 1) {
    #     # plot single gene as barplot
    p1 <-
        gene_exp %>%
        mutate(facet = str_c('Top ', top_n)) %>% 
        slice(seq.int(n(), by = -1L, length.out = top_n)) %>% 
        ggplot(aes(x=tissue, y=expression)) + 
        geom_col(aes(fill=tissue)) +
        geom_text(aes(y = ymx*0.04, label = tissue), hjust=0) +
        labs(title = str_c("GTEx: ", gene)) +
        theme_bw() + 
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.y = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.x = element_blank(),
              plot.margin = margin(t = 2)) + 
        guides(fill='none') +
        ylim(0, ymx) +
        coord_flip() +
        facet_wrap(~facet, strip.position = 'left')
    
    p2 <-
        gene_exp %>%
        mutate(facet = 'Other Tissues') %>% 
        slice(seq.int(n()-top_n)) %>% 
        ggplot(aes(x=tissue, y=expression)) + 
        geom_col(aes(), col = 'gray35', fill = 'gray34') +
        labs(y = "log-2 TPM") +
        theme_bw() + 
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.y = element_blank(),
              plot.margin = margin()) + 
        ylim(0, ymx) +
        coord_flip() +
        facet_wrap(~facet, strip.position = 'left')
    
    
    p <- plot_grid(p1, p2, rel_heights = c(1,2), ncol = 1,
                   align = 'v')
    
    return(p)
}





