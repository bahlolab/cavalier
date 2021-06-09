
# creates IGV snapshots using xvfb-run and igv.sh
# Can either set igv_sh to full path to igv.sh, or Sys.setenv(PATH=...), likewise for xvfb-run/ singularity
#' @importFrom purrr walk map_chr
#' @importFrom stringr str_c
#' @export
create_igv_snapshots <- function(candidates, bams, genome,
                                 output_dir = 'igv_snapshots',
                                 overwrite = FALSE,
                                 slop = 20,
                                 igv_sh = 'igv.sh',
                                 width = 720,
                                 height = 720,
                                 xvfb_run = 'xvfb-run',
                                 singularity_img = NULL,
                                 singularity_bin = 'singularity',
                                 genome_file = NULL,
                                 name_panel_width = 10,
                                 prefs = character()) {

    snapshot_tbl <-
        candidates %>% 
        transmute(chrom = chromosome,
                  start = position,
                  end = start + nchar(reference)) %>% 
        mutate(start = start - slop,
               end = end + slop,
               filename = file.path(output_dir, str_c(chrom, start, end, sep='_')) %>% 
                   str_c('.png'))
    candidates$igv_filename <- snapshot_tbl$filename
    
    if (!overwrite) {
        snapshot_tbl <- filter(snapshot_tbl, !file.exists(filename))
    }
    
    if (nrow(snapshot_tbl)) {
        # create output dir
        if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }
        # create igv user dir and genome dir
        user_dir <- file.path(output_dir, 'igv')
        genome_dir <- file.path(user_dir, 'genomes')
        dir.create(genome_dir, showWarnings = FALSE, recursive = TRUE)
        # link genome
        if (!is.null(genome_file)) {
            if (!file.exists(genome_link)){
                file.symlink(normalizePath(genome_file), genome_dir)
            }
        }
        # write prefs.properties to set IGV dimensions
        c(str_c('NAME_PANEL_WIDTH=', name_panel_width),
          str_c('IGV.Bounds=0,0,', width, ',', height),
          prefs) %>%
            write_lines(file.path(user_dir, 'prefs.properties'))
        # write batch script
        batchfile <- file.path(output_dir, 'igv_snapshots.batch')
        c('new',
          str_c('genome ', genome),
          map_chr(bams, ~ str_c('load ', .)),
          str_c('maxPanelHeight ', height - 95),
          pmap(snapshot_tbl, function(chrom, start, end, filename, ...) {
              c(str_c('goto ', chrom, ':', start, '-', end),
                str_c('snapshot ', filename))
          }) %>% unlist(),
          'exit') %>% 
            write_lines(batchfile)
        
        # snapshot command
        cmd <-
            str_c(xvfb_run,
                  '--auto-servernum',
                  '--server-num=1',
                  str_c('-s \'-screen 0 ', width, 'x', height, 'x8\''),
                  igv_sh,
                  '--igvDirectory', user_dir,
                  '-b', batchfile,
                  sep = ' ')
        # wrap singularity if using
        if (!is.null(singularity_img)) {
            bind_dirs <- 
                c(getwd(),
                  output_dir,
                  normalizePath(bams) %>% dirname()) %>% 
                normalizePath() %>% 
                unique() %>% 
                remove_child_dirs()
            cmd <-
                str_c(singularity_bin, 
                      'exec',
                      str_c('-B ', bind_dirs, collapse = ' '),
                      singularity_img,
                      cmd,
                      sep = ' ')
        }
        message('executing: ', cmd)
        # execute command
        stopifnot(system(cmd) == 0)
    }
    
    return(candidates)
}
