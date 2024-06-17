
#' Get assembly gap and centromere locations for get_cavalier_opt('ref_genome')
#' 
#' @importFrom dplyr select mutate across bind_rows arrange_all filter
#' @importFrom readr read_tsv cols
get_centromeres_gaps <- function() {
  
  ref_genome <- get_cavalier_opt('ref_genome')
  
  fun <- function() {
    
    if (ref_genome == 'hg38') {
      agp_url <- get_cavalier_opt("agp_url_hg38")
    } else {
      agp_url <- get_cavalier_opt("agp_url_hg19")
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
      
      cen_url <- get_cavalier_opt("cen_url_hg38")
      
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
  }
  
  cache(
    fun = fun,
    name = str_c('UCSC.assembly_gaps.', ref_genome)
  )
}

