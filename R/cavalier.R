#' cavalier: Candidate Variant List Evaluation Report
#'
#' \pkg{cavalier} assists and automates variant interpretation in next-generation sequencing data:
#' 1) filter list of variants and create a short list of candidate variants
#' 2) combine information from multiple sources useful for evaluating and prioritising variants
#' 3) create PDF output for visualing and communicating results
#'
#' @docType package
#' @name cavalier

NULL


# environment to store default options, user settable with function cavalier_options()
cavalier_opts <- new.env()
# set default options
cavalier_opts$cache_dir <- '~/.cavalier'
cavalier_opts$snapshot_dir <- '.igv_snapshots'
cavalier_opts$retry_pause_base <- 5
cavalier_opts$retry_pause_min <- 5
cavalier_opts$retry_times <- 5
cavalier_opts$ref_genome <- 'hg38'
cavalier_opts$xvfb_run_cmd <- 'xvfb-run'
cavalier_opts$igv_cmd <- 'igv.sh'
cavalier_opts$singularity_cmd <- 'singularity'
cavalier_opts$igv_hg38_uri <- 'https://s3.amazonaws.com/igv.org.genomes/hg38/hg38.genome'
cavalier_opts$igv_hg19_uri <- 'https://s3.amazonaws.com/igv.org.genomes/hg19/hg19.genome'
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
  walk2(names(dots), dots, function(n, v) {
    assign(n, v, envir = cavalier_opts)
  })
} 


