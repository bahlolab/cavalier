
#' @importFrom stringr str_c str_extract str_split
#' @importFrom magrittr '%>%' set_colnames
#' @importFrom dplyr bind_cols
#' @importFrom tidyr unnest chop
#' @importFrom tidyr unnest chop
#' @importFrom SeqArray seqSetFilter
load_vcf <- function(vcf_filename, 
                     samples = NULL,
                     info_columns = c('AF', 'AC', 'DP', 'QD', 'MQ', 'FS', 'SOR', 'MQRankSum', 'ReadPosRankSum', 'InbreedingCoeff')) 
{
  gds <- vcf_to_gds(vcf_filename)
  
  # restrict to biallelic variants
  num_allele <- SeqArray::seqNumAllele(gds)
  if (!all(num_allele == 2L)) {
    warning('Multiallelic sites present - these will be ignored. ',
            'Please first flatten VCF using "bcftools norm -m-any --do-not-normalize <VCF>" or equivalent. ',
            'Using ', sum(num_allele == 2, na.rm = T), ' biallelic sites of ', length(num_allele), ' total sites.')
    seqSetFilter(gds, num_allele == 2)
  }
  
  vcf_data <-
    bind_cols(get_vcf_fixed(gds),
              get_vcf_info(gds, info_columns) %>% select(-variant_id)) %>% 
    mutate(genotype = get_vcf_genotypes(gds, samples)) %>% 
    bind_cols(get_vcf_sample_AD(gds, samples)) %>% 
    mutate(GQ = get_vcf_sample_GQ(gds, samples))
  
  return(vcf_data)
}

#' @importFrom SeqArray seqOpen seqVCF2GDS
#' @importFrom assertthat assert_that
#' @importFrom rlang is_scalar_character
#' @importFrom stringr str_replace
vcf_to_gds <- function(vcf_fn) {
  
  assert_that(is_scalar_character(vcf_fn), 
             file.exists(vcf_fn))
  
  gds_fn <- str_replace(vcf_fn, '.vcf.gz', '.gds')
  if (! file.exists(gds_fn)) {
    seqVCF2GDS(vcf.fn = vcf_fn, out.fn = gds_fn, storage.option = 'ZIP_RA', ignore.chr.prefix = '')
  }
  seqOpen(gds_fn, allow.duplicate = T)
}

#' @importFrom dplyr tibble mutate
#' @importFrom SeqArray seqGetData
get_vcf_fixed <- function(gds) {
  tibble(
    variant_id = seqGetData(gds, 'variant.id'),
    chrom = seqGetData(gds, 'chromosome'),
    pos = seqGetData(gds, 'position'),
    qual = SeqArray::qual(gds),
    filt = SeqArray::filt(gds)) %>% 
    bind_cols(seqGetData(gds, 'allele') %>% 
                str_split_fixed(',',2) %>% 
                set_colnames(c('ref', 'alt')) %>% 
                as_tibble())
}

#' @importFrom dplyr tibble mutate as_tibble
#' @importFrom SeqArray seqGetData
get_vcf_info <- function(gds, info_columns) {
  tibble(
    variant_id = seqGetData(gds, 'variant.id')) %>% 
    bind_cols(str_c('annotation/info/', info_columns) %>% 
                setNames(info_columns) %>% 
                map_dfc(function(x) {
                  seqGetData(gds, x) %>% 
                    { `if`(is.list(.), .$data, .) }
                })
    )
}

#' @importFrom dplyr everything
#' @importFrom magrittr set_colnames set_rownames
get_vcf_genotypes <- function(gds, samples) {
  
  if (!is.null(samples)) {
    seqSetFilter(gds, sample.id = samples)
  }
  
  genotypes <- seqGetData(gds, 'genotype')
  str_c(pmin(genotypes[1, ,], genotypes[2, ,]),
        pmax(genotypes[1, ,], genotypes[2, ,]),
        sep = '/') %>% 
    matrix(nrow = dim(genotypes)[2], ncol = dim(genotypes)[3]) %>% 
    t() %>% 
    set_colnames(seqGetData(gds, 'sample.id')) %>% 
    as_tibble()
}

#' @importFrom dplyr as_tibble
#' @importFrom magrittr set_colnames set_rownames
get_vcf_sample_AD <- function(gds, samples) {
  
  if (!is.null(samples)) {
    seqSetFilter(gds, sample.id = samples)
  }
  
  AD <- seqGetData(gds, 'annotation/format/AD')
  assert_that(all(AD$length == 2))
  nvar <- ncol(AD$data) / 2
  
  setNames(1:2, c('depth_ref', 'depth_alt')) %>% 
    map(function(i) {
      (AD$data[, seq.int(from = i, by = 2, length.out = nvar)]) %>% 
        t() %>% 
        set_colnames(seqGetData(gds, 'sample.id')) %>%
        as_tibble()
    }) %>% 
    as_tibble()
}

#' @importFrom dplyr as_tibble
#' @importFrom magrittr set_colnames set_rownames
get_vcf_sample_GQ <- function(gds, samples) {
  
  if (!is.null(samples)) {
    seqSetFilter(gds, sample.id = samples)
  }
  
  GQ <- seqGetData(gds, 'annotation/format/GQ')
  assert_that(all(GQ$length == 1))

  t(GQ$data) %>% 
    set_colnames(seqGetData(gds, 'sample.id')) %>% 
    as_tibble()
}
