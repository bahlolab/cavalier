#' @export
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
annotate_gaps <- function(variants,
                          min_overlap_bp = 100,
                          min_overlap_prop = 0.25) {
  
  gap_ranges <-
    get_centromeres_gaps() %>% 
    na.omit() %>% 
    with(GRanges(chrom, IRanges(start, end), type = type))
  
  var_ranges <-
    variants %>% 
    mutate(row = row_number(),
           chrom = chrom,
           start = pos,
           end = if_else(is.na(END), pos, END)) %>%
    select(row, chrom, start, end) %>% 
    na.omit() %>% 
    filter(start <= end) %>% 
    with(GRanges(chrom, IRanges(start, end), row = row))
  
  matches <-
    suppressWarnings(
    GenomicRanges::findOverlaps(var_ranges, gap_ranges, minoverlap = min_overlap_bp) %>% 
    as_tibble() %>% 
    chop(-queryHits)  %>% 
      mutate(overlap = map2_dbl(queryHits, subjectHits, function(vh, gh) {
        sum(
          GenomicRanges::width(
            GenomicRanges::intersect(var_ranges[vh], gap_ranges[gh]))) / 
          GenomicRanges::width(var_ranges[vh])
      }))) %>% 
      filter(overlap >= min_overlap_prop) %>% 
    unchop(subjectHits) %>% 
    transmute(row = var_ranges[queryHits]$row,
              gap_type = gap_ranges[subjectHits]$type) %>% 
    chop(gap_type) %>%
    mutate(gap_type = map_chr(gap_type, function(gt) {
      sort(gt) %>% unique() %>% str_c(collapse = '&')
    }))
  
  variants %>% 
    mutate(row = row_number()) %>% 
    left_join(matches, by = 'row') %>% 
    select(-row)
}