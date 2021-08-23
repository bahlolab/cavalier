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

# environment to store cached tables, e.g. gtex expression, omim genemap etc
cavalier_cache <- new.env()

# environment to store default options, user settable with function cavalier_options()
cavalier_opts <- new.env()
# set default options
cavalier_opts$cache_dir <- '~/.cavalier'
cavalier_opts$igv_snapshot_dir <- '.igv_snapshots'
cavalier_opts$ref_genome <- 'hg38'
cavalier_opts$xvfb_run_cmd <- 'xvfb-run'
cavalier_opts$igv_cmd <- 'igv.sh'
cavalier_opts$singularity_cmd <- 'singularity'
cavalier_opts$hgnc_complete_uri <- 'http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/monthly/tsv/hgnc_complete_set_2021-08-01.txt'
cavalier_opts$rvis_uri <- 'http://genic-intolerance.org/data/RVIS_Unpublished_ExACv2_March2017.txt'
cavalier_opts$gevir_uri <- 'http://www.gevirank.org/static/files/gene_ranking.csv'
cavalier_opts$gtex_gene_median_tpm_uri <- 'https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz'
cavalier_opts$igv_hg38_uri <- 'https://s3.amazonaws.com/igv.org.genomes/hg38/hg38.genome'
cavalier_opts$igv_hg19_uri <- 'https://s3.amazonaws.com/igv.org.genomes/hg19/hg19.genome'
cavalier_opts$gtex_tissues <-
  c("Brain - Amygdala",
    "Brain - Anterior cingulate cortex (BA24)",
    "Brain - Caudate (basal ganglia)",
    "Brain - Cerebellar Hemisphere",
    "Brain - Cerebellum",
    "Brain - Cortex",
    "Brain - Frontal Cortex (BA9)",
    "Brain - Hippocampus",
    "Brain - Hypothalamus",
    "Brain - Nucleus accumbens (basal ganglia)",
    "Brain - Putamen (basal ganglia)",
    "Brain - Spinal cord (cervical c-1)",
    "Brain - Substantia nigra",
    "Adipose - Subcutaneous",
    "Artery - Tibial",
    "Breast - Mammary Tissue",
    "Esophagus - Mucosa",
    "Lung",
    "Muscle - Skeletal",
    "Nerve - Tibial",
    "Skin - Not Sun Exposed (Suprapubic)",
    "Skin - Sun Exposed (Lower leg)",
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


