#' Return OMIM phenotype and inhertance model for given gene symbol
#' 
#' @param gene gene symbol
#' @param genemap2 location of genemap2.txt file downloaded from OMIM (see https://omim.org/downloads/)
#' @return data frame of OMIM phenotypes and inheritance associated with gene
# #' @examples
# #' ***TODO***

omim_table <- function(gene, genemap2=NULL)
{
    if ((is.null(genemap2) | !file.exists(genemap2)) & !("genemap2_omim_table" %in% ls())) {
        print("Warning: no OMIM genemap2 table found.")
        print("Download genemap2.txt from https://omim.org/downloads/ or specify location of file.")
        return(NULL)
    } else if (!("genemap2_omim_table" %in% ls())) {
        require(tidyverse)
        genemap2_omim_table <- 
            read.delim(genemap2, skip=3, stringsAsFactors=FALSE) %>% 
            mutate(alt_symbol = map(Gene.Symbols, 
                                 ~ str_trim(c(str_split(., ',', simplify = TRUE))))) %>% 
            unnest(alt_symbol) %>% 
            as.data.frame()
        assign("genemap2_omim_table", genemap2_omim_table, envir=.GlobalEnv)
    }

    # Find OMIM matches for gene
    gene_gm2 <- genemap2_omim_table[genemap2_omim_table$alt_symbol %in% gene, ]
    if (nrow(gene_gm2) == 0) {
        return(NULL)
    }
        
    phenotypes <- paste(gene_gm2$Phenotypes, collapse="; ")
    if (phenotypes == "") {
        return(NULL)
    }

    pt_split <- strsplit(phenotypes, "; ")[[1]]
    pt_split2 <- strsplit(pt_split, ")")
    pt_pheno <- sapply(pt_split2, function(x){paste0(paste0(ifelse(length(x) == 1, x, x[-length(x)]), collapse=")"), ")")})
    pt_pheno <- gsub("?", "", gsub("(2)", "", gsub("(3)", "", pt_pheno, fixed=TRUE), fixed=TRUE), fixed=TRUE)
    pt_inher <- sapply(pt_split2, function(x){ifelse(length(x) > 1, x[length(x)], "NA")})
    # *** BELOW IF STATEMENTS NEED TO BE TESTED MORE THOROUGHLY ***
    if (length(pt_inher) > 0) {
        pt_inher <- ifelse(startsWith(pt_inher, ", "), substr(pt_inher, 3, 999), pt_inher)
    }
    if (is.character(pt_inher)) {
        pt_inher <- gsub(", ", "\n", pt_inher, fixed=TRUE)
    }
    OMIMtable <- data.frame(OMIM.phenotype=pt_pheno, OMIM.inheritance=pt_inher, stringsAsFactors=FALSE)
    if (nrow(OMIMtable) == 0 | ncol(OMIMtable) == 0) {
        return(NULL)
    }
    OMIMtable$OMIM.inheritance[OMIMtable$OMIM.inheritance == "NA"] <- ""
    rownames(OMIMtable) <- NULL
    colnames(OMIMtable) <- c("OMIM phenotype", "inheritance")
    if (all(OMIMtable$inheritance == "")) {
        OMIMtable <- OMIMtable[, "OMIM phenotype", drop=FALSE]
        OMIMtable[, "OMIM phenotype"] <- sapply(OMIMtable[, "OMIM phenotype"], function(x){paste0(strwrap(x, width=60), collapse="  \n")})
    } else {
        OMIMtable[, "OMIM phenotype"] <- sapply(OMIMtable[, "OMIM phenotype"], function(x){paste0(strwrap(x, width=42), collapse="  \n")})
    }
    return(OMIMtable)
}
