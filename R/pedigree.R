
#' @importFrom assertthat assert_that
#' @importFrom rlang is_scalar_character
#' @importFrom readr read_tsv
#' @export
read_ped <- function(file) 
{
  ped_df <- 
    read_tsv(file, 
             col_names = c('fid', 'iid', 'pid', 'mid', 'sex', 'phe'),
             col_types = 'ccccii') %>% 
    mutate(pid = if_else(pid == '0', NA_character_, pid),
           mid = if_else(mid == '0', NA_character_, mid))
  
  ped_df %>% 
    with(assert_that(all(pid %in% iid | is.na(pid)),
                     all(mid %in% iid | is.na(mid)),
                     length(iid) == length(unique(iid)),
                     all(sex %in% c('1', '2')),
                     all(phe %in% c('1', '2'))))
  
  ped_df 
}

# identify trios in families
# only returns an unambgiuous trio, e.g. a single affected child with two parents per family
#' @importFrom dplyr filter select group_by ungroup n
#' @export
get_trio <- function(ped_df) 
{
  assert_that(is.data.frame(ped_df))
  ped_df %>% 
    filter(phe == 2, !is.na(pid), !is.na(mid)) %>% 
    select(fid, iid, pid, mid, phe, sex) %>% 
    group_by(fid) %>% 
    filter(n() == 1) %>% 
    ungroup()
}

#' @export
get_affected <- function(ped_df) 
{
  filter(ped_df, phe ==2) %>% pull(iid)
}

#' @export
get_unaffected <- function(ped_df) 
{
  filter(ped_df, phe == 1) %>% pull(iid)
}

#' @export
plot_ped <- function(ped_df,
                     draw = TRUE,
                     cex = 1.3,
                     col = 'darkred',
                     mar = c(0,2,1,2)) 
{
  pedigree <- as_kinship_pedigree(ped_df)
  
  lab <- `if`('label' %in% colnames(ped_df),
              ped_df$label, ped_df$iid)
  
  grob <- cowplot::as_grob(
    ~ kinship2::plot.pedigree(pedigree, lab,
                              mar = mar,
                              col = col,
                              cex = cex)
  )
  
  if (draw) {
    grid::grid.newpage()
    grid::grid.draw(grob)
    invisible(grob)
  } else {
    grob
  }
}

# note: assumes single family as plotting seems broken with multiple families
as_kinship_pedigree <- function(ped_df)
{
  ped_df %>% 
    with(
      kinship2::pedigree(id = iid,
                         dadid = pid, 
                         momid = mid, 
                         sex = sex,
                         affected = phe))
}


