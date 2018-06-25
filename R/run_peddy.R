#' Run peddy python tool for fast pedigree, sex and ancestry check
#' 
#' Pedersen & Quinlan, 2017, The American Journal of Human Genetics, 100, 3, 406-413
#' Who’s Who? Detecting and Resolving Sample Anomalies in Human DNA Sequencing Studies with Peddy
#' https://doi.org/10.1016/j.ajhg.2017.01.017
#' https://github.com/brentp/peddy

#' @param output_dir cavalier output directory
#' @param pedigree pedigree (.ped) file location
#' @param vcf VCF file location
#' @param processors number of processors to use when running peddy (default: 1)
# #' @examples
# #' ***TODO***

run_peddy <- function(output_dir, pedigree, vcf, processors=1)
{
    if (!dir.exists(paste0(output_dir, "data/peddy/"))) {
    	dir.create(paste0(output_dir, "data/peddy/"), showWarnings=TRUE, recursive=TRUE)
    }
    command <- paste("python -m peddy --plot -p", processors, "--prefix peddy", vcf, pedigree)
    curr_wd <- getwd()
    setwd(paste0(endslash_dirname(output_dir), "data/peddy/"))
    system(command)
    setwd(curr_wd)
}
