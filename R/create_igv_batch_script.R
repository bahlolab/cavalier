#' Create IGV batch script
#' 
#' Creates igv_batch.txt in [output_dir]/data/ directory, which when run in IGV
#' after manually loading BAM files creates snapshots in [output_dir]/data/igv_output/
#' 
#' @param candidates candidate variants data.frame
#' @param output_dir output base directory
#' @param slop number of base pairs to include on either side of variant (default: 10)
#' @return candidate variants data.frame with additional 'igv_filename' column containing filename of IGV snapshot
# #' @examples
# #' ***TODO***

create_igv_batch_script <- function(candidates, output_dir, slop=10)
{
    output_dir <- endslash_dirname(output_dir)
    data_dir <- paste0(output_dir, "data/")

    if (!dir.exists(data_dir)) {
        dir.create(data_dir, showWarnings=TRUE, recursive=TRUE)
    }

    igv_script_filename <- paste0(data_dir, "igv_batch.txt")
    igv_output_dir <- paste0(data_dir, "igv_output/")

    # Create IGV filename for each candidate variant
    candidates$igv_filename <- paste0(igv_output_dir, candidates$chromosome, "_", candidates$start, 
                                    ifelse(candidates$start == candidates$end, "", paste0("_", candidates$end)),
                                    ".png")

    igv_script <- paste("snapshotDirectory", igv_output_dir, "\n")
    for (ii in 1:nrow(candidates)) {
        ii_chr <- candidates[ii, "chromosome"]
        ii_start <- candidates[ii, "start"]
        ii_end <- candidates[ii, "end"]
        ii_filename <- basename(as.character(candidates[ii, "igv_filename"]))

        igv_script <- paste0(igv_script, "goto ", ii_chr, ":", ii_start - slop, "-", ii_end + slop, "\n")
        igv_script <- paste0(igv_script, "snapshot ", ii_filename, "\n")
    }

    cat(igv_script, file=igv_script_filename)
    return(candidates)
}

