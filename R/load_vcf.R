
#' @importFrom stringr str_c str_extract str_split
#' @importFrom magrittr '%>%' set_colnames
#' @importFrom dplyr bind_cols left_join
#' @importFrom tidyr unnest chop
#' @importFrom tidyr unnest chop
#' @importFrom SeqArray seqSetFilter
#' @importFrom rlang is_bool is_scalar_character
#' @export
load_vcf <- function(input, 
                     samples = NULL,
                     info_columns = c('AF', 'AC', 'QD'),
                     annot_source = 'VEP',
                     annot_source_field = 'CSQ',
                     add_annot = default_annotations(),
                     remove_chrom_prefix = NULL,
                     add_chrom_prefix = NULL) 
{
  # check args
  assert_that(is_gds(input) | (is_scalar_character(input) && file.exists(input)),
              is.null(samples) | is.character(samples),
              is.null(info_columns) | is.character(info_columns),
              is.null(annot_source) | is_scalar_character(annot_source),
              is.null(annot_source) | is_scalar_character(annot_source_field),
              is.character(add_annot) & all(add_annot %in% all_annotations()),
              is.null(remove_chrom_prefix) | is_scalar_character(remove_chrom_prefix),
              is.null(add_chrom_prefix) | is_scalar_character(add_chrom_prefix),
              is.null(remove_chrom_prefix) | is.null(add_chrom_prefix))
  
  gds <- `if`(is_gds(input), input, vcf_to_gds(input))

  # restrict to biallelic variants
  num_allele <- SeqArray::seqNumAllele(gds)
  if (!all(num_allele == 2L)) {
    warning('Multiallelic sites present - these will be ignored. ',
            'Please first flatten VCF using "bcftools norm -m-any --do-not-normalize <VCF>" or equivalent. ',
            'Using ', sum(num_allele == 2, na.rm = T), ' biallelic sites of ', length(num_allele), ' total sites.')
    seqSetFilter(gds, num_allele == 2)
  }
  
  if (!is.null(samples)) {
    seqSetFilter(gds, sample.id = samples)
  }
  
  variants <-
    bind_cols(get_vcf_fixed(gds),
              get_vcf_info(gds, info_columns) %>% select(-variant_id)) %>% 
    mutate(genotype = get_vcf_genotypes(gds)) %>% 
    bind_cols(get_vcf_sample_AD(gds)) %>% 
    mutate(GQ = get_vcf_sample_GQ(gds)) %>% 
    mutate(chrom = case_when(
      !is.null(add_chrom_prefix)    ~ str_c(add_chrom_prefix, str_remove(chrom, str_c('^', add_chrom_prefix))),
      !is.null(remove_chrom_prefix) ~ str_remove(chrom, str_c('^', remove_chrom_prefix)),
      TRUE                          ~ chrom)) %>% 
    (function(data) {
      `if`(annot_source == 'VEP',
           data %>% 
             left_join(get_vep_ann(gds, annot_source_field, add_annot = add_annot),
                       by = 'variant_id'),
           data)
    })
  
  return(variants)
}

#' @export
default_annotations <- function() {
  c('rvis_percentile', 'gevir_percentile', 'grantham_score')
}

#' @export
all_annotations <- function() {
  c('rvis_percentile', 'gevir_percentile', 'loeuf_percentile', 'grantham_score')
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
  seqOpen(gds_fn, allow.duplicate = TRUE)
}

is_gds <- function(object) {
  inherits(object, 'SeqVarGDSClass')
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
#' @importFrom purrr map_dfc

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
get_vcf_genotypes <- function(gds) {
  
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
get_vcf_sample_AD <- function(gds) {
  
  AD <- seqGetData(gds, 'annotation/format/AD')
  assert_that(all(AD$length == 2))
  nvar <- ncol(AD$data) / 2
  
  setNames(1:2, c('depth_ref', 'depth_alt')) %>% 
    map(function(i) {
      (AD$data[, seq.int(from = i, by = 2, length.out = nvar), drop = FALSE]) %>% 
        t() %>% 
        set_colnames(seqGetData(gds, 'sample.id')) %>%
        as_tibble()
    }) %>% 
    as_tibble()
}

#' @importFrom dplyr as_tibble
#' @importFrom magrittr set_colnames set_rownames
get_vcf_sample_GQ <- function(gds, samples) {
  
  GQ <- seqGetData(gds, 'annotation/format/GQ')
  assert_that(all(GQ$length == 1))

  t(GQ$data) %>% 
    set_colnames(seqGetData(gds, 'sample.id')) %>% 
    as_tibble()
}
