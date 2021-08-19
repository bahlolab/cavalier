
#' @importFrom assertthat assert_that
#' @importFrom rlang is_scalar_character is_scalar_double
#' @importFrom purrr reduce
#' @importFrom magrittr and
#' @export
add_inheritance <- function(variants,
                            af_column,
                            affected = NULL,
                            unaffected = NULL,
                            pedigree = NULL,
                            models = inheritance_models(),
                            af_dominant = 1e-4,
                            af_recessive = 1e-2,
                            af_compound_het = 1e-2,
                            af_de_novo = 1e-4)
{
  assert_that(is.data.frame(variants),
              is.data.frame(variants$genotype),
              is_scalar_character(af_column),
              is.character(affected) | is.null(affected),
              is.character(unaffected) | is.null(unaffected),
              !is.null(affected) & !is.null(pedigree),
              is.null(affected) & is.null(pedigree),
              is.character(models) & length(models) > 0,
              is_scalar_character(pedigree) | is.null(pedigree),
              is_scalar_double(af_dominant),
              is_scalar_double(af_recessive),
              is_scalar_double(af_compound_het),
              is_scalar_double(af_de_novo))
  
  
  genotype <- as.matrix(variants$genotype)
  if (is.null(pedigree)) {
    
    if (is.null(unaffected)) { 
      unaffected <- character() 
    }
    assert_that(all(affected %in% colnames(genotype)),
                all(unaffected %in% colnames(genotype)))
    
  } else {
    stop('pedigree not yet implemented')
  }
  
  
  inh_df <- select(variants, variant_id)
  # dominant
  if ('dominant' %in% models) {
    inh_df <- inh_df %>% 
      mutate(dominant = 
               c(map(affected,   ~ genotype[[.]] %in% c('0/1', '1/1')),
                 map(unaffected, ~ genotype[[.]] == c('0/0')),
                 list(replace_na(variants[[af_column]] <= af_dominant, TRUE))) %>% 
               reduce(and))
  }
  # recessive
  if ('recessive' %in% models) {
    inh_df <- inh_df %>% 
      mutate(recessive = 
               c(map(affected,   ~ genotype[[.]] == '1/1'),
                 map(unaffected, ~ genotype[[.]] != c('1/1')),
                 list(replace_na(variants[[af_column]] <= af_recessive, TRUE))) %>% 
               reduce(and))
  }
  # compound_het
  if ('compound_het' %in% models) {
    inh_df <- inh_df %>% 
      mutate(recessive = 
               c(map(affected,   ~ genotype[[.]] == '1/1'),
                 map(unaffected, ~ genotype[[.]] != c('1/1')),
                 list(replace_na(variants[[af_column]] <= af_recessive, TRUE))) %>% 
               reduce(and))
  }
  # de novo
  if ('de_novo' %in% models) {
    warning('de_novo inheritance model not yet implemented')
  }
}


inheritance_models <- function() {
  c('dominant', 'recessive', 'compound_het', 'de_novo')
}