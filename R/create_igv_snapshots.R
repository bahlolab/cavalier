#' Create IGV snapshots
#' 
#' Creates IGV images automatically using IGV-snapshot-automator
#' 
#' @param candidates candidate variants data.frame
#' @param bams bam file(s) to load into IGV and take snapshots of
#' @param reference_genome reference genome to use (hg19, hg38)
#' @param output_dir output base directory
#' @param slop number of base pairs to include on either side of variant (default: 20)
#' @param overwrite overwrite existing files with same name (default: FALSE)
#' @return candidate variants data.frame with additional 'igv_filename' column containing filename of IGV snapshot
# #' @examples
# #' ***TODO***

# *** TOFIX: location of IGV-snapshot-automator is hardcoded to lab_bahlo share ***

make_IGV_snapshots <- "/wehisan/bioinf/lab_bahlo/software/apps/IGV-snapshot-automator/IGV-snapshot-automator_2019-02-21/make_IGV_snapshots.py"
IGV_jar <- "/wehisan/bioinf/lab_bahlo/software/apps/IGV-snapshot-automator/IGV-snapshot-automator_2019-02-21/bin/IGV_2.3.81/igv.jar"


create_igv_snapshots <- function(candidates, bams, reference_genome, output_dir, slop=20, overwrite=FALSE)
{
    output_dir <- endslash_dirname(output_dir)
    igv_output_dir <- paste0(output_dir, "data/igv_output/")
    
    if (!dir.exists(igv_output_dir)) {
        dir.create(igv_output_dir, showWarnings=TRUE, recursive=TRUE)
    }
    
    # Create IGV snapshot filenames and bed file for IGV-snapshot automator
    length_ref <- sapply(candidates$reference, nchar)
    length_alt <- sapply(candidates$alternate, nchar)
    start <- candidates$position
    end <- ifelse(length_ref <= length_alt, start, start + length_ref - length_alt - 1)
    
    candidates$igv_filename <- paste0(igv_output_dir, candidates$chromosome, "_", start - slop, "_", end + slop, "_h500.png")
    
    bed <- data.frame(candidates$chromosome, start - slop, end + slop, candidates$igv_filename, stringsAsFactors=FALSE)
    if (!overwrite) {
        bed <- bed[!file.exists(candidates$igv_filename), ]
    }
    
    if (nrow(bed) > 0) {
        temp_igv_bed <- paste0(igv_output_dir, "temp_igv.bed")
        write.table(bed, file=temp_igv_bed, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
        
        command <- paste0("python ", make_IGV_snapshots, 
                            " ", paste(bams, collapse=" "), 
                            " -g ", reference_genome, 
                            " -bin ", IGV_jar, 
                            " -r ", temp_igv_bed,
                            " -o ", igv_output_dir)
        system(command)
        file.remove(temp_igv_bed)
    }
    return(candidates)
}


