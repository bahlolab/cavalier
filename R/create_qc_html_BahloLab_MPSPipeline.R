#' Create QC HTML summary reports from Bahlo Lab inhouse pipeline output
#' 
#' @param MPSPipeline_output_dir location of Bahlo lab MPSPipeline output directory
#' @param output_dir output base directory
# #' @examples
# #' ***TODO***

# *** TODO: generalise this function for other pipelines
create_qc_html_BahloLab_MPSPipeline <- function(MPSPipeline_output_dir, output_dir)
{
    MPSPipeline_output_dir <- endslash_dirname(MPSPipeline_output_dir)
    output_dir <- endslash_dirname(output_dir)

    ### Create FastQC summary table

    # Copy FastQC dir from location in MPSPipeline to html_files
    system(paste0("cp -r ", MPSPipeline_output_dir, "fastq_summary/FastQC/ ", output_dir, "html_files/fastqc/"))

    fastqc_summary <- read.csv(paste0(output_dir, "html_files/fastqc/fastQCmultiple_sample_summary/summary_table.csv"), stringsAsFactors=FALSE)
    
    summary_cols <- c("BSQ", "TSQ", "SQS", "BSC", "SGCC", "BNC", "SLD", "OS", "AC", "KMER")

    fastqc_summary <- fastqc_summary[, c("Filename", summary_cols)]
    
    fastqc_summary_colours <- fastqc_summary
    fastqc_summary_colours[, 1] <- "black"
    for (cc in summary_cols) {
        fastqc_summary_colours[, cc][fastqc_summary[, cc] == "PASS"] <- "green"
        fastqc_summary_colours[, cc][fastqc_summary[, cc] == "WARN"] <- "orange"
        fastqc_summary_colours[, cc][fastqc_summary[, cc] == "FAIL"] <- "red"
    }

    # Replace text with specific link for each cell
    col_imagenames <- c("", "per_base_quality", "per_tile_quality", "per_sequence_quality", "per_base_sequence_content", "per_sequence_gc_content", "per_base_n_content", "sequence_length_distribution", "", "adapter_content", "")
    for (ii in 1:nrow(fastqc_summary)) {
        ii_imagedir <- paste0("fastqc/", gsub(".fastq", "_fastqc", gsub(".fastq.gz", "_fastqc", as.character(fastqc_summary[ii, 1]), fixed=TRUE), fixed=TRUE), "/Images/")
        for (jj in 2:ncol(fastqc_summary)) {
            jj_imagename <- col_imagenames[jj]
            if (jj_imagename != "") {
                fastqc_summary[ii, jj] <- paste0("<a href='", ii_imagedir, jj_imagename, ".png'> <span style='color:", fastqc_summary_colours[ii, jj],"'>", fastqc_summary[ii, jj], "</span></a>")
            } else {
                fastqc_summary[ii, jj] <- paste0("<span style='color:", fastqc_summary_colours[ii, jj], "'>", fastqc_summary[ii, jj], "</span>")
            }
        }
    }
    fastqc_summary$Report <- paste0("<a href='fastqc/", gsub(".fastq", "_fastqc", gsub(".fastq.gz", "_fastqc", as.character(fastqc_summary[, 1]), fixed=TRUE), fixed=TRUE), ".html'>report</a>")

    # Output FastQC table
    dt <- DT::datatable(fastqc_summary, escape=FALSE)
    DT::saveWidget(dt, paste0(output_dir, "html_files/fastqc_summary.html"), selfcontained=FALSE)


    ### Create alignment stats table

    align_file <- list.files(path=MPSPipeline_output_dir, pattern="*summarystats*", full.names=TRUE)
    align_file <- align_file[!endsWith(align_file, ".stderr.txt")]
    align_stats <- read.delim(align_file, sep="\t")
    colnames(align_stats) <- c(colnames(align_stats)[-1], "none")
    align_summary <- data.frame(sampleID=unique(sapply(rownames(align_stats), function(x)strsplit(x, "_")[[1]][1])), stringsAsFactors=FALSE)
    rownames(align_summary) <- align_summary$sampleID
    for (sID in rownames(align_summary)) {
        if ("Unique" %in% colnames(align_stats)) {
            sID_table <- align_stats[startsWith(rownames(align_stats), sID), ]
            align_summary[sID, "Total Reads"] <- sum(sID_table$All)
            align_summary[sID, "Uniquely Aligned"] <- sID_table$Unique[1]
            align_summary[sID, "UA %"] <- round(100 * align_summary[sID, "Uniquely Aligned"] / align_summary[sID, "Total Reads"], 1)
            align_summary[sID, "Non-duplicate uniquely aligned"] <- sID_table$Duprem[1]
            align_summary[sID, "ND UA %"] <- round(100 * align_summary[sID, "Non-duplicate uniquely aligned"] / align_summary[sID, "Total Reads"], 1)
        } else if ("UniqueAlignment" %in% colnames(align_stats)) {
            sID_table <- align_stats[startsWith(rownames(align_stats), sID), ]
            align_summary[sID, "Total Reads"] <- sID_table$ReadSequences[1]
            align_summary[sID, "Uniquely Aligned"] <- sID_table$UniqueAlignment[1]
            align_summary[sID, "UA %"] <- round(100 * align_summary[sID, "Uniquely Aligned"] / align_summary[sID, "Total Reads"], 1)
            align_summary[sID, "Non-duplicate uniquely aligned"] <- sID_table$Duprem[1]
            align_summary[sID, "ND UA %"] <- round(100 * align_summary[sID, "Non-duplicate uniquely aligned"] / align_summary[sID, "Total Reads"], 1)
        }
    }

    # Output alignment summary table
    dt <- DT::datatable(align_summary, escape=FALSE)
    DT::saveWidget(dt, paste0(output_dir, "html_files/align_stats.html"), selfcontained=FALSE)


    ### Create coverage stats table  (*** TODO: add coverage plots ***)

    coverage_file <- list.files(path=paste0(MPSPipeline_output_dir, "coverage/capture"), pattern="*sample_summary*", full.names=TRUE)
    coverage_stats <- read.delim(coverage_file, sep="\t")
    coverage_stats <- coverage_stats[-c(nrow(coverage_stats)), c(1, 6, 5, 3, 4, 8, 9)]
    colnames(coverage_stats) <- c("sampleID", "1st quartile", "median", "mean", "3rd quartile", "% bases > 5", "% bases > 10")

    # Output coverage summary table
    dt <- DT::datatable(coverage_stats, escape=FALSE)
    DT::saveWidget(dt, paste0(output_dir, "html_files/coverage_stats.html"), selfcontained=FALSE)
}

