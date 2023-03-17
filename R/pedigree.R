
#' @importFrom assertthat assert_that
#' @importFrom rlang is_scalar_character
#' @importFrom readr read_tsv
#' @export
read_ped <- function(file) {
  
  ped_df <-
    read_tsv(
      file, 
      col_names = c('famid', 'id', 'dadid', 'momid', 'sex', 'affected'),
      col_types = 'ccccii') %>% 
    distinct()
  
  message('read pedigree with ',
          nrow(ped_df),
          ' individuals and ',
          n_distinct(ped_df$famid),
          ' families')
  
  ped_df <-
    ped_df %>% 
    mutate(
      across(c(dadid, momid), ~ if_else(. == '0', NA_character_, .)),
      sex = case_when(
        sex == 1 ~ 'male',
        sex == 2 ~ 'female',
        TRUE     ~ 'unknown'),
      affected = case_when(
        affected == 2 ~ TRUE,
        affected == 1 ~ FALSE,
        TRUE          ~ NA)
    ) %>% 
    ped_fix_sex() %>% 
    ped_add_parents()
  
  # check validity
  ped_df %>% 
    with(assert_that(
      all(dadid %in% id | is.na(dadid)),
      all(momid %in% id | is.na(momid)),
      length(id) == length(unique(id)),
      all(sex %in% c('male', 'female', 'unknown')),
      all(is.logical(affected))))
  
  ped_df
}

#' @export
#' @importFrom readr write_tsv
write_ped <- function(...) 
{
  write_tsv(..., col_names = F)
}


ped_fix_sex <- function(ped_df) {
  # check no contradictions
  mom_and_dad <- with(ped_df, momid[!is.na(momid) & momid %in% dadid])
  if (length(mom_and_dad)) {
    rlang::abort(
      str_c("individual is recored as both mother and father: ",
            str_c(mom_and_dad, collapse = ', ')))
  }
  # make sure all dads are male and moms are female
  male_moms <-
    ped_df %>% 
    filter(sex == 'male',
           id %in% momid) %>% 
    pull(id)
  
  female_dads <-
    ped_df %>% 
    filter(sex == 'female',
           id %in% dadid) %>% 
    pull(id)
  
  if (length(male_moms)) {
    message('converting to female as recorded as mother: ', 
            str_c(male_moms, collapse = ', '))
  }
  
  if (length(female_dads)) {
    message('converting to male as recorded as father: ', 
            str_c(female_dads, collapse = ', '))
  }
  
  ped_df %>% 
    mutate(sex = case_when(
      id %in% female_dads ~ 'male',
      id %in% male_moms   ~ 'female',
      TRUE                ~ sex
    ))
}

ped_add_parents <- function(ped_df,
                            mom_suff = '-mother',
                            dad_suff = '-father') {
  
  miss_dad <- with(ped_df, which(is.na(dadid) & !is.na(momid)))
  miss_mom <- with(ped_df, which(is.na(momid) & !is.na(dadid)))
  
  res <-
    ped_df %>% 
    mutate(
      dadid = replace(dadid, miss_dad, str_c(id[miss_dad], dad_suff)),
      momid = replace(momid, miss_mom, str_c(id[miss_mom], mom_suff)))
  
  miss_dad <- with(res, setdiff(dadid, id)) %>% na.omit()
  miss_mom <- with(res, setdiff(momid, id)) %>% na.omit()
  
  missing <-
    inner_join(
      bind_rows(
        select(res, famid, id = dadid),
        select(res, famid, id = momid)) %>% 
        filter(id %in% c(miss_dad, miss_mom)) %>% 
        distinct(),
      tibble(id = c(miss_dad, miss_mom),
             sex = c(rep('male', length(miss_dad)), rep('female', length(miss_mom))),
             dadid = NA_character_,
             momid = NA_character_,
             affected = NA
      ),
      by = 'id'
    )
  
  if (nrow(missing)) {
    message('added ', nrow(missing), ' missing parental records')
  }
  
  bind_rows(res, missing) %>% 
    arrange_all()
}



# identify trios in families
# only returns an unambgiuous trio, e.g. a single affected child with two parents per family
#' @importFrom dplyr filter select group_by ungroup n
#' @export
get_trio <- function(ped_df) 
{
  assert_that(is.data.frame(ped_df))
  ped_df %>% 
    filter(affected, !is.na(dadid), !is.na(momid)) %>% 
    select(famid, id, dadid, momid, affected, sex) %>% 
    group_by(famid) %>% 
    filter(n() == 1) %>% 
    ungroup()
}

#' @export
get_affected <- function(ped_df) 
{
  filter(ped_df, affected) %>% pull(id)
}

#' @export
get_unaffected <- function(ped_df) 
{
  filter(ped_df, !affected) %>% pull(id)
}

#' @importFrom cowplot as_grob
#' @importFrom ggplotify as.ggplot
#' @export
plot_ped <- function(ped_df,
                     cex = 1.3,
                     col = 'darkred',
                     mar = c(0,2,1,2)) 
{
  pedigree <- as_kinship_pedigree(ped_df)
  
  lab <- `if`('label' %in% colnames(ped_df),
              ped_df$label, ped_df$id)
  
  grob <- cowplot::as_grob(
    ~ kinship2::plot.pedigree(pedigree, lab,
                              mar = mar,
                              col = col,
                              cex = cex)
  )
  as.ggplot(grob)
}

as_kinship_pedigree <- function(ped_df, single_family = TRUE) 
{
  if (single_family) {
    with(ped_df,
         kinship2::pedigree(id = id,
                            dadid = dadid, 
                            momid = momid, 
                            sex = sex,
                            affected = affected))
  } else {
    with(ped_df,
         kinship2::pedigree(id = id,
                            dadid = dadid, 
                            momid = momid, 
                            sex = sex,
                            affected = affected,
                            famid = famid))
  }
}

