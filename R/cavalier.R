#' cavalier: Candidate Variant List Evaluation Report
#'
#' \pkg{cavalier} assists and automates variant interpretation in next-generation sequencing data:
#' 1) filter list of variants and create a short list of candidate variants
#' 2) combine information from multiple sources useful for evaluating and prioritising variants
#' 3) create PDF output for visualing and communicating results
#'
#' @docType package
#' @name cavalier
#' @export HGNC_alias

NULL

options('cavalier.cache_dir' = '~/.cavalier')
options('cavalier.hgnc_complete_uri' = 'http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/monthly/tsv/hgnc_complete_set_2021-08-01.txt')
options('cavalier.gtex_gene_median_tpm_uri' = 'https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz')
options('cavalier.rvis_uri' = 'http://genic-intolerance.org/data/RVIS_Unpublished_ExACv2_March2017.txt')


clear_cache <- function()
{
  cache_dir <- getOption('cavalier.cache_dir')
  if (dir.exists(cache_dir)) {
    file.remove(list.files(cache_dir, full.names = TRUE))
  }
}

cache <- function(fun, fn) 
{
  if (file.exists(fn)) {
    readRDS(fn)
  } else {
    res <- fun()
    tmp_fn <- tempfile(pattern = basename(fn) %>% str_c('.'),
                       tmpdir = dirname(fn))
    saveRDS(res, tmp_fn)
    file.rename(tmp_fn, fn)
    res
  }
}
