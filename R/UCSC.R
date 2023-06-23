
#' @importFrom dplyr select mutate across bind_rows arrange_all filter
#' @importFrom readr read_tsv cols

get_centromeres_gaps <- function() {
  
  ref_genome <- get_cavalier_opt('ref_genome')
  
  (function() {

    if (ref_genome == 'hg38') {
      agp_url <- 'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.agp.gz'
    } else {
      agp_url <- 'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.agp.gz'
    }
    
    gaps <- 
      read_tsv(agp_url,
               col_names = c('object',
                             'object_beg',
                             'object_end',
                             'part_number',
                             'component_type',
                             'gap_length',
                             'gap_type',
                             'linkage',
                             'linkage_evidence'),
               col_types = cols()) %>% 
      filter(component_type %in% c('N', 'U')) %>% 
      select(chrom = object,
             start = object_beg,
             end = object_end,
             type = gap_type)  %>% 
      arrange_all()
    
    if (ref_genome == 'hg38') {
      
      cen_url <- 'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz'
      
      gaps <- 
        bind_rows(
          gaps,
          read_tsv(cen_url,
                   col_names = c('NUM',
                                 'chrom',
                                 'start',
                                 'end',
                                 'id'),
                   col_types = cols()) %>% 
            mutate(type = 'centromere') %>% 
            select(chrom, start, end, type)
        ) %>% 
        arrange_all()
    }
    
    return(gaps)
  }) %>% 
    cache(str_c(ref_genome, '.centromeres_gaps'), disk = TRUE)
}

