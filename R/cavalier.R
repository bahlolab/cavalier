#' cavalier: Candidate Variant List Evaluation Report
#'
#' \pkg{cavalier} assists and automates variant interpretation in next-generation sequencing data:
#' 1) filter list of variants and create a short list of candidate variants
#' 2) combine information from multiple sources useful for evaluating and prioritising variants
#' 3) create PDF output for visualing and communicating results
#'
#' @name cavalier
"_PACKAGE"


# environment to store default options, user settable with function cavalier_options()
cavalier_opts <- new.env()
# set default options
cavalier_opts$cache_dir <- '~/.cavalier'
cavalier_opts$snapshot_dir <- '.igv_snapshots'
cavalier_opts$retry_pause_base <- 5
cavalier_opts$retry_pause_min <- 5
cavalier_opts$retry_times <- 3
cavalier_opts$ref_genome <- 'hg38'
cavalier_opts$xvfb_run_cmd <- 'xvfb-run'
cavalier_opts$igv_cmd <- 'igv.sh'
cavalier_opts$singularity_cmd <- 'singularity'
cavalier_opts$igv_hg38_uri <- 'https://s3.amazonaws.com/igv.org.genomes/hg38/hg38.genome'
cavalier_opts$igv_hg19_uri <- 'https://s3.amazonaws.com/igv.org.genomes/hg19/hg19.genome'
############ GTEX options #################
# note: gtex_gene_median_tpm_url can be set to local file with e.g. "file:///path/to/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
cavalier_opts$gtex_gene_median_tpm_url <-
  "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
cavalier_opts$gtex_tissues <-
  c(
    "Adipose - Subcutaneous",
    "Artery - Tibial",
    "Bladder",
    "Breast - Mammary Tissue",
    "Brain - Cerebellum",
    "Brain - Cortex",
    "Colon - Sigmoid",
    "Esophagus - Mucosa",
    "Heart - Left Ventricle",
    "Kidney - Cortex",
    "Liver",
    "Lung",
    "Muscle - Skeletal",
    "Nerve - Tibial",
    "Ovary",
    "Stomach",
    "Skin - Sun Exposed (Lower leg)",
    "Testis",
    "Thyroid",
    "Whole Blood")

############ HGNC options #################
# either "latest", "local" or monthly release e.g. "2024-06-04"
cavalier_opts$hgnc_ver = "latest"
cavalier_opts$hgnc_monthly_base_url <- 'http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/monthly/tsv/'
# set to a local file path for use without web access
cavalier_opts$hgnc_local_file <- NULL

#' @export
get_cavalier_opt <- function(name = NULL) {
  if (is.null(name)) {
    return(as.list(cavalier_opts))
  }
  cavalier_opts[[name]]
}
 
#' @export
set_cavalier_opt <- function(...) {
  dots <- dots_list(...)
  assert_that(is_named(dots))
  purrr::walk2(names(dots), dots, function(n, v) {
    assign(n, v, envir = cavalier_opts)
  })
}
