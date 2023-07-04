
#' @importFrom assertthat assert_that
#' @importFrom rlang is_scalar_character is_scalar_double is_scalar_integerish
#' @importFrom purrr reduce map_chr
#' @importFrom magrittr and not
#' @importFrom dplyr select mutate filter summarise_all across all_of summarise first
#' @importFrom tidyr pivot_longer complete chop unchop
#' @importFrom stringr str_c
#' @export
add_inheritance <- function(variants,
                            af_column = 'af_gnomad',
                            uid_column = 'hgvs_genomic',
                            ped_file = NULL,
                            affected = NULL,
                            unaffected = NULL,
                            models = inheritance_models(),
                            min_depth = 5L,
                            af_dominant = 1e-4,
                            af_recessive = 1e-2,
                            af_compound_het = 1e-2,
                            af_de_novo = 1e-4,
                            compound_het_max = 3)
{
  assert_that(is.data.frame(variants),
              is.data.frame(variants$genotype),
              is.null(min_depth) || is.data.frame(variants$depth_ref),
              is.null(min_depth) || is.data.frame(variants$depth_alt),
              is_scalar_character(af_column),
              all(c('variant_id', af_column) %in% colnames(variants)),
              is.character(affected) | is.null(affected),
              is.character(unaffected) | is.null(unaffected),
              xor(is.null(ped_file), is.null(affected)),
              is.character(models) & length(models) > 0,
              is_scalar_character(ped_file) | is.null(ped_file),
              is.null(min_depth) || is_scalar_integerish(min_depth),
              is_scalar_double(af_dominant),
              is_scalar_double(af_recessive),
              is_scalar_double(af_compound_het),
              is_scalar_double(af_de_novo))
  
  
  genotype <- variants$genotype
  sample_set <- colnames(variants$genotype)
  
  if (is.null(ped_file)) {
    
    if (is.null(unaffected)) { 
      unaffected <- character() 
    }
    trio <- tibble()
    
  } else {
    
    ped_df <- read_ped(ped_file)
    affected <- get_affected(ped_df)
    n_aff <- length(affected)
    affected <- intersect(affected, sample_set)
    unaffected <- get_unaffected(ped_df) %>% intersect(sample_set)
    trio <- get_trio(ped_df) %>% 
      filter(id %in% sample_set,
             dadid %in% sample_set,
             momid %in% sample_set)
  }
  
  assert_that(all(affected %in% sample_set),
              all(unaffected %in% sample_set),
              length(affected) > 0)
  
  sample_set <- c(affected, unaffected)
  
  if (!is.null(min_depth)) {
    sample_min_depth <- 
      (variants$depth_ref[,sample_set] + 
         variants$depth_alt[,sample_set]) %>% 
      apply(1, min)
    gte_min_depth <- sample_min_depth >= min_depth
  } else {
    gte_min_depth <- rep(TRUE, nrow(variants))
  }
  
  
  inh_df <- tibble(id = seq_len(nrow(variants)), pair_uid = NA_character_)
  # dominant
  if ('dominant' %in% models) {
    inh_df <- inh_df %>% 
      mutate(dominant = 
               c(map(affected,   ~ genotype[[.]] %in% c('0/1', '1/1')),
                 map(unaffected, ~ genotype[[.]] == c('0/0'))) %>% 
               reduce(and) %>% 
               and(replace_na(variants[[af_column]] <= af_dominant, TRUE)) %>% 
               and(gte_min_depth))
  }
  # recessive
  if ('recessive' %in% models) {
    inh_df <- inh_df %>% 
      mutate(recessive = 
               c(map(affected,   ~ genotype[[.]] == '1/1'),
                 map(unaffected, ~ genotype[[.]] != c('1/1'))) %>% 
               reduce(and) %>% 
               and(replace_na(variants[[af_column]] <= af_recessive, TRUE)) %>% 
               and(gte_min_depth))
  }
  # compound_het
  if ('compound_het' %in% models) {
    gene_gt <- 
      tibble(id = inh_df$id,
             gene = variants$gene) %>% 
      bind_cols(genotype[, sample_set]) %>% 
      filter(gte_min_depth,
             replace_na(variants[[af_column]] <= af_compound_het, TRUE)) %>% 
      na.omit()
    
    c_h_df <-
      gene_gt %>% 
      select(id, gene, all_of(affected)) %>% 
      pivot_longer(c(-gene, -id), 
                   names_to = 'sample',
                   values_to = 'gt') %>% 
      filter(gt == '0/1') %>% 
      select(-gt) %>% 
      chop(id) %>% 
      filter(lengths(id) > 1,
             lengths(id) <= compound_het_max) %>% 
      (function(data) {
        # split cases of more than two sites into multiples
        bind_rows(filter(data, lengths(id) == 2) %>% 
                    mutate(id = as.list(id)),
                  filter(data, lengths(id) > 2) %>% 
                    mutate(id = map(id, function(id) {
                      combn(id, 2) %>% as.data.frame() %>% as.list() %>% unname()
                    })) %>% 
                    unnest(id))
      }) %>% 
      (function(data) {
        # ensure other samples also consistent
        `if`(length(sample_set) > 1,
             filter(data, map2_lgl(sample, id, function(sample, vid) {
               gene_gt %>% 
                 filter(id %in% vid) %>% 
                 select(all_of(setdiff(sample_set, sample))) %>% 
                 summarise_all(~ all(. == '0/1') | any(. == '1/1')) %>% 
                 mutate(across(all_of(unaffected), not)) %>% 
                 unlist() %>% 
                 all()
             })),
             data)
      }) %>% 
      select(id) %>% 
      distinct() %>% 
      mutate(
             pair_id = id) %>% 
      unchop(id) %>% 
      mutate(pair_id = map2_int(pair_id, id, setdiff)) %>% 
      left_join(select(variants, pair_uid = all_of(uid_column)) %>% 
                  mutate(pair_id = inh_df$id),
                by = 'pair_id') %>% 
      group_by(id) %>% 
      summarise(compound_het = TRUE,
                pair_uid = str_c(pair_uid, collapse = '&'),
                .groups = 'drop')
    
    inh_df <-
      left_join(inh_df, c_h_df, by = c('id', 'pair_uid')) %>% 
      mutate(compound_het = replace_na(compound_het, FALSE))
  }
  # de novo
  if ('de_novo' %in% models && nrow(trio) && length(affected) == 1) {
    inh_df <- inh_df %>% 
      mutate(de_novo = 
               c(map(affected,   ~ genotype[[.]] == '0/1'),
                 map(unaffected, ~ genotype[[.]] == c('0/0'))) %>% 
               reduce(and) %>% 
               and(replace_na(variants[[af_column]] <= af_de_novo, TRUE)) %>% 
               and(gte_min_depth))
  }
  
  inh_df <-
    inh_df %>% 
    pivot_longer(where(is.logical), 
                 names_to = 'inheritance') %>% 
    filter(value) %>% 
    select(-value) %>% 
    chop(-id) %>% 
    (function(data) {
      `if`('compound_het' %in% models,
           mutate(data,
                  pair_uid = map_chr(pair_uid,
                                     ~ `if`(all(is.na(.)), NA_character_,
                                            first(na.omit(.))))),
           data)
    }) %>% 
    complete(id = seq_len(nrow(variants)),
             fill = list(inheritance = vctrs::list_of(NA_character_))) %>% 
    mutate(inheritance = map_chr(inheritance, ~ str_c(., collapse = '&'))) %>% 
    arrange(id)
  
  bind_cols(variants,
            select(inh_df, -id))
}

inheritance_models <- function() {
  c('dominant', 'recessive', 'compound_het', 'de_novo')
}
