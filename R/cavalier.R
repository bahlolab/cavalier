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
cavalier_opts$retry_pause_base <- 1
cavalier_opts$retry_pause_cap  <- 30
cavalier_opts$retry_times <- 3
cavalier_opts$ref_genome <- 'hg38'
cavalier_opts$xvfb_run_cmd <- 'xvfb-run'
cavalier_opts$igv_cmd <- 'igv.sh'
cavalier_opts$singularity_cmd <- 'singularity'
cavalier_opts$igv_hg38_uri <- 'https://s3.amazonaws.com/igv.org.genomes/hg38/hg38.genome'
cavalier_opts$igv_hg19_uri <- 'https://s3.amazonaws.com/igv.org.genomes/hg19/hg19.genome'
cavalier_opts$use_memoisation <- TRUE

# three modes: 
#   "latest" - use latest version from web, fail if can't access it
#   "fallback" - use latest version if from web if possible, fallback to latest cache version if not
#   "offline" - use latest cached version, only check online if no cached version available
cavalier_opts$database_mode <- "fallback"

########### Gene Intolerance #############
# gene intolerance urls, TODO: add gnomadv4 intolerances
cavalier_opts$gevir_url <- "http://www.gevirank.org/static/files/gene_ranking.csv" 

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
# NULL or specific monthly release e.g. "2024-06-04" or "local" for local file
cavalier_opts$hgnc_ver = NULL
cavalier_opts$hgnc_monthly_base_url <- 'http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/monthly/tsv/'
# set to a local file path for use without web access
cavalier_opts$hgnc_local_file <- NULL

########### HPO options ####################
cavalier_opts$hpo_api_base_url <- "https://ontology.jax.org/api/"
cavalier_opts$hpo_api_max_failuers <- 10L
cavalier_opts$hpo_github_url <- "https://github.com/obophenotype/human-phenotype-ontology/"

########## PanelApp options ###############
cavalier_opts$panelapp_urls <- list(
  PAA = "https://panelapp.agha.umccr.org/",       # PanelApp Australia
  PAE = "https://panelapp.genomicsengland.co.uk/" # PanelApp England
)

#' @export
get_cavalier_opt <- function(name = NULL) {
  if (is.null(name)) {
    return(as.list(cavalier_opts))
  }
  cavalier_opts[[name]]
}

######### OMIM options ####################
cavalier_opts$mim2gene_url <- 'https://omim.org/static/omim/data/mim2gene.txt'

######## UCSC options #####################
cavalier_opts$agp_url_hg38 <- 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.agp.gz'
cavalier_opts$cen_url_hg38 <- 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz'
cavalier_opts$agp_url_hg19 <- 'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.agp.gz'

 
#' @export
set_cavalier_opt <- function(...) {
  dots <- dots_list(...)
  assert_that(is_named(dots))
  purrr::walk2(names(dots), dots, function(n, v) {
    assign(n, v, envir = cavalier_opts)
  })
}
