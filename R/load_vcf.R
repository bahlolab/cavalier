
#' @importFrom stringr str_c str_extract str_split
#' @importFrom magrittr '%>%' set_colnames
#' @importFrom dplyr bind_cols left_join
#' @importFrom tidyr unnest chop
#' @importFrom tidyr unnest chop
#' @importFrom SeqArray seqSetFilter seqClose
#' @importFrom rlang is_bool is_scalar_character
#' @export
load_vcf <- function(input, 
                     samples = NULL,
                     caller = 'GATK',
                     annotater = 'VEP',
                     annotater_field = 'CSQ',
                     info_columns = caller_info_columns(caller),
                     format_columns = caller_format_columns(caller),
                     additional_annotation = default_annotations(),
                     SVO = FALSE,
                     remove_chrom_prefix = NULL,
                     add_chrom_prefix = NULL) 
{
  # check args
  assert_that(is_gds(input) | (is_scalar_character(input) && file.exists(input)),
              is.null(samples) | is.character(samples),
              is_scalar_character(caller) && is_valid_caller(caller),
              is.null(annotater) | is_scalar_character(annotater) && is_valid_annotater(annotater),
              is_scalar_character(annotater_field),
              is.character(info_columns),
              is.character(format_columns),
              is.character(additional_annotation),
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
    get_vcf_fixed(gds) %>% 
    bind_cols(get_vcf_info(gds, info_columns) %>%
                select(-variant_id)) %>% 
    (function(x) {
      `if`(length(seqGetData(gds, 'sample.id')),
           mutate(x, genotype = get_vcf_genotypes(gds)),
           x)
    }) %>%  
    (function(x) {
      `if`('AD' %in% format_columns,
           bind_cols(x, get_vcf_sample_AD(gds)),
           x) 
    }) %>% 
    (function(x) {
      `if`('GQ' %in% format_columns,
           mutate(x, GQ = get_vcf_sample_GQ(gds)),
           x) 
    }) %>% 
    mutate(chrom = case_when(
      !is.null(add_chrom_prefix)    ~ str_c(add_chrom_prefix, str_remove(chrom, str_c('^', add_chrom_prefix))),
      !is.null(remove_chrom_prefix) ~ str_remove(chrom, str_c('^', remove_chrom_prefix)),
      TRUE                          ~ chrom)) %>% 
    (function(x) {
      `if`(annotater == 'VEP',
           left_join(x, get_vep_ann(gds, annotater_field, add_annot = additional_annotation, SVO = SVO),
                     by = 'variant_id'),
           x)
    })
  
  seqClose(gds)
  
  return(variants)
}

is_valid_caller <- function(x) { 
  x %in% c('GATK', 'Mutect2', 'manta', 'manta-jasmine') 
}

is_valid_annotater <- function(x) { 
  x %in% c('VEP', '') 
}

caller_info_columns <- function(caller) {
  
  assert_that(is_scalar_character(caller))
  
  if (caller == 'GATK') {
    c('AF', 'AC', 'AN', 'QD')
  } else if (caller == 'Mutect2') {
    c('AF', 'AC', 'AN')
  } else if (caller == 'manta') {
    c('AF', 'AC', 'AN', 'END', 'SVTYPE', 'SVLEN')
  } else if (caller == 'manta-jasmine') {
    c('AF', 'AC', 'AN', 'END', 'SVTYPE', 'SVLEN', 'AVG_LEN', 'AVG_START', 'AVG_END')
  }
}

caller_format_columns <- function(caller) {
  
  assert_that(is_scalar_character(caller))
  
  if (caller == 'GATK') {
    c('AD', 'GQ')
  } else if (caller == 'Mutect2') {
    c('AD')
  } else if (caller == 'manta') {
    character()
  } else if (caller == 'manta-jasmine') {
    character()
  }
}

#' @export
default_annotations <- function() {
  c('loeuf_percentile', 'gevir_percentile', 'grantham_score')
}

#' @export
all_annotations <- function() {
  c('gevir_percentile', 'loeuf_percentile', 'grantham_score')
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
    seqVCF2GDS(vcf.fn = vcf_fn,
               out.fn = gds_fn,
               storage.option = 'ZIP_RA',
               ignore.chr.prefix = '',
               verbose = FALSE)
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
    id = seqGetData(gds, 'annotation/id'),
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
  
  var_id <- seqGetData(gds, 'variant.id')
  cols <- 
    setNames(info_columns, info_columns) %>% 
    map_dfc(function(icol) {
      
      x <- 
        tryCatch(seqGetData(gds, str_c('annotation/info/', icol)),
                 error = function(e) { NULL })
      if (is.null(x)) {
        warning('INFO/', icol, ' does not exists, returning NA')
        return(NA)
      }
      
      if (is.list(x)) {
        if (all(x$length) == 1) {
          val <- x$data
        } else {
          val <- vector(mode(x$data), length(x$length))
          first <- setdiff(unique(cumsum(x$length)), 0)
          val[x$length == 0] <- as(NA, mode(x$data))
          val[x$length != 0] <- x$data[first]
          if (any(x$length) > 1) {
            warning('INFO/', icol, ' has entries with length > 1, using first value')
          }
        }
      } else {
        val <- x
      }
      if (length(val) != length(var_id)) {
        warning('INFO/', icol, ' had length ', length(val), ' but expected ', length(var_id), ', returning NA') 
        val <- NA
      }
      return(val)
    })
  
  bind_cols(variant_id = var_id, cols) %>% 
    (function(x) {
      cols <- names(x) %>% intersect(c('AVG_LEN', 'AVG_START', 'AVG_END'))
      `if`(length(cols),
           mutate(x, across(all_of(cols), ~ as.integer(round(as.numeric(.))))),
           x)
    })
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
  
  # depending on SeqArray version, sometimes get a list and sometimes a matrix
  if (is.list(GQ)) {
    assert_that(all(GQ$length == 1))
    GQ <- GQ$data
  }
  
  t(GQ) %>% 
    set_colnames(seqGetData(gds, 'sample.id')) %>% 
    as_tibble()
}
