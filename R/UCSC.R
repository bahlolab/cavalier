
#' @importFrom dplyr select mutate across bind_rows arrange_all filter

get_centromeres_gaps <- function() {
  
  ref_genome <- get_cavalier_opt('ref_genome')
  
  (function() {
    mySession <- rtracklayer::browserSession()
    rtracklayer::genome(mySession) <- ref_genome
    
    gaps <- 
      rtracklayer::getTable(rtracklayer::ucscTableQuery(mySession, track="gap")) %>% 
      as_tibble() %>% 
      select(chrom, start = chromStart, end=chromEnd, type) %>% 
      mutate(across(where(is.factor), as.character),
             start = start + 1L) %>% 
      mutate(type = if_else(type %in% c('centromere', 'heterochromatin', 'short_arm', 'telomere'),
                            type, 'gap'))
    
    if (ref_genome == 'hg38') {
      centro <- 
        rtracklayer::getTable(rtracklayer::ucscTableQuery(mySession, track="centromeres")) %>% 
        as_tibble() %>% 
        select(chrom, start = chromStart, end=chromEnd) %>% 
        mutate(type = 'centromere') %>% 
        mutate(across(where(is.factor), as.character),
               start = start + 1L)
      
      gaps <-
        bind_rows(gaps, centro) %>% 
        arrange_all()
    } else {
      gaps <-
        gaps %>% 
        filter(str_detect(chrom, '^chr[0-9XY]+$')) %>% 
        arrange_all()
    }
    
    return(gaps)
  }) %>% 
    cache(str_c(ref_genome, '.centromeres_gaps'), disk = TRUE)
}
