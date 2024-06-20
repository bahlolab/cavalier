
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
                            gene_column = 'ensembl_gene_id',
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
                            compound_het_max = Inf)
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
  
  
  inh_df <- tibble(
    id = seq_len(nrow(variants)),
    uid = variants[[uid_column]],
    gene = variants[[gene_column]],
  )
  
  # dominant
  if ('dominant' %in% models) {
    inh_df <- 
      inh_df %>% 
      mutate(dominant = 
               c(map(affected,   ~ genotype[[.]] %in% c('0/1', '1/1')),
                 map(unaffected, ~ genotype[[.]] == c('0/0'))) %>% 
               reduce(and) %>% 
               and(replace_na(variants[[af_column]] <= af_dominant, TRUE)) %>% 
               and(gte_min_depth))
  }
  # recessive
  if ('recessive' %in% models) {
    inh_df <- 
      inh_df %>% 
      mutate(recessive = 
               c(map(affected,   ~ genotype[[.]] == '1/1'),
                 map(unaffected, ~ genotype[[.]] != c('1/1'))) %>% 
               reduce(and) %>% 
               and(replace_na(variants[[af_column]] <= af_recessive, TRUE)) %>% 
               and(gte_min_depth))
  }
  # compound_het
  if ('compound_het' %in% models) {
    
    get_ch_cand <- function(data) {
      data %>% 
        na.omit() %>% 
        distinct() %>% 
        pivot_longer(-c(uid, gene), 
                     names_to = 'sample',
                     values_to = 'gt') %>% 
        mutate(n_allele = str_count(gt, '1')) %>% 
        select(-gt) %>% 
        distinct() %>% 
        chop(-c(gene, sample)) %>% 
        mutate(n_allele_gene = map_int(n_allele, sum)) %>% 
        group_by(gene) %>% 
        filter(all(n_allele_gene >= 2)) %>% 
        ungroup() %>% 
        select(gene, uid, sample) %>% 
        mutate(uid = map(uid, function(x) {
          # split into pairs
          combn(x, 2) %>% as.data.frame() %>% as.list() %>% unname()
        })) %>% 
        unnest(uid)
    }
    # need > 2 alleles per sample in a gene
    cand_aff <-
      inh_df %>% 
      select(uid, gene, recessive) %>% 
      bind_cols(genotype[, affected]) %>% 
      filter(replace_na(variants[[af_column]] <= af_compound_het, TRUE)) %>% 
      filter(!recessive) %>% 
      select(-recessive) %>% 
      get_ch_cand()
    
    if (length(unaffected)) {
      # find any present in unaffected individuals
      cand_unaff <-
        inh_df %>% 
        select(uid, gene) %>% 
        bind_cols(genotype[, unaffected]) %>% 
        filter(uid %in% unlist(cand_aff$uid)) %>% 
        get_ch_cand() %>% 
        select(gene, uid) %>% 
        distinct()
      # remove from cand_aff
      cand_aff <-
        cand_aff %>% 
        anti_join(cand_unaff) %>% 
        group_by(gene) %>% 
        filter(n_distinct(sample) == length(affected))
    }
    
    cand_aff <-
      cand_aff %>% 
      unnest(uid) %>% 
      select(gene, uid) %>% 
      distinct() %>% 
      mutate(compound_het = TRUE) %>% 
      group_by(gene) %>% 
      mutate(pair_uid = map_chr(uid, ~ str_c(setdiff(uid, .), collapse = '&'))) %>% 
      ungroup()
      
    
    inh_df <-
      inh_df %>% 
      left_join(
        cand_aff,
        by = c('uid', 'gene')
      ) %>% 
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
    select(-uid, -gene) %>% 
    pivot_longer(where(is.logical), 
                 names_to = 'inheritance') %>% 
    bind_rows(tibble(pair_uid = character())) %>% 
    filter(value) %>% 
    select(-value) %>% 
    chop(-id)  %>% 
    complete(
      id = seq_len(nrow(variants)),
      fill = list(inheritance = vctrs::list_of(NA_character_),
                  pair_uid = vctrs::list_of(NA_character_))) %>% 
    mutate(inheritance = map_chr(inheritance, ~ str_c(., collapse = '&')),
           pair_uid = map_chr(pair_uid, first)) %>% 
    arrange(id)
  
  bind_cols(variants, select(inh_df, inheritance, pair_uid))
}

inheritance_models <- function() {
  c('dominant', 'recessive', 'compound_het', 'de_novo')
}
