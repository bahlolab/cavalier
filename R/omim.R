
#' @importFrom tidyr replace_na separate_rows
#' @importFrom stringr str_c str_remove str_remove_all str_extract
#' @importFrom dplyr "%>%" mutate rename select if_else summarise group_by bind_rows arrange_all pull first
#' @export
get_omim_table <- function(genemap2)
{
    if (!file.exists(genemap2)) {
        stop("Provided OMIM genemap2 table does not exist. ", 
             "Download genemap2.txt from https://omim.org/downloads/ and specify file location.")
    } 
    
    omim_table <- getOption('cavalier.omim_table')
    
    if (is.null(omim_table)) {
        
        col_names <- c(
            'chromosome', 'genomic_position_start', 'genomic_position_end', 'cyto_location', 
            'computed_cyto_location', 'mim_number', 'gene_symbols', 'gene_name', 
            'approved_gene_symbol', 'entrez_gene_id', 'ensembl_gene_id', 'comments', 
            'phenotypes', 'mouse_gene_symbol_id')
        
        hgnc_sym <- get_hgnc_complete() %>% pull(symbol)
        
        omim_table <- 
            read_tsv(genemap2, col_names = col_names, comment = '#', col_types = cols()) %>% 
            select(gene_symbols, ensembl_gene_id, entrez_gene_id, phenotypes) %>% 
            mutate(entrez_gene_id = as.integer(entrez_gene_id)) %>% 
            filter(!is.na(phenotypes)) %>% 
            # match symbol with ensemble_gene_id and entrez_gene_id first
            mutate(symbol = hgnc_ensembl2sym(ensembl_gene_id),
                   symbol = if_else(is.na(symbol), hgnc_entrez2sym(entrez_gene_id), symbol)) %>% 
            # next try hgnc_symbol
            (function(data) {
                filter(data, is.na(symbol)) %>% 
                    mutate(id = seq_along(gene_symbols)) %>% 
                    separate_rows(gene_symbols, sep = ',\\s+') %>% 
                    mutate(symbol_1 = if_else(gene_symbols %in% hgnc_sym, gene_symbols, NA_character_),
                           symbol_2 = hgnc_sym2sym(gene_symbols, remove_unknown = TRUE)) %>% 
                    group_by(id, phenotypes) %>% 
                    summarise(symbol = if_else(any(!is.na(symbol_1)),
                                               na.omit(symbol_1) %>% first(),
                                               na.omit(symbol_2) %>% first()),
                    .groups = 'drop') %>% 
                    filter(!is.na(symbol)) %>% 
                    select(-id) %>% 
                    bind_rows(select(data, symbol, ensembl_gene_id, entrez_gene_id, phenotypes), .) %>% 
                    arrange_all()
            }) %>% 
            separate_rows(phenotypes, sep=';\\s+') %>% 
            mutate(phenotype = str_extract(phenotypes, '^.*(?= \\([0-9]\\))') %>% 
                       str_remove_all('[\\[\\]?{}]'),
                   inheritance = str_extract(phenotypes, '(?<= \\([0-9]\\), ).*$')) %>%
            select(-phenotypes)
        
        
        options('cavalier.omim_table' = omim_table)
    }
    
    return(omim_table)
}
