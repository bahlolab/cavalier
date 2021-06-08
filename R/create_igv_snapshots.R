
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
                                 igv_max_height = 650,
                                 igv_args = character(),
                                 xvfb_run = 'xvfb-run',
                                 singularity_img = NULL,
                                 singularity_bin = 'singularity',
                                 width = 720,
                                 height = 500) {

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
        # write batch script
        if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }
        batch_tmp <- tempfile(tmpdir = output_dir, pattern = 'igv_', fileext = '.bat')
        c('new',
          str_c('genome ', genome),
          map_chr(bams, ~ str_c('load ', .)),
          str_c('maxPanelHeight ', igv_max_height),
          pmap(snapshot_tbl, function(chrom, start, end, filename, ...) {
              c(str_c('goto ', chrom, ':', start, '-', end),
                str_c('snapshot ', filename))
          }) %>% unlist(),
          'exit') %>% 
            write_lines(batch_tmp)
        
        # snapshot command
        cmd <-
            str_c(xvfb_run,
                  '--auto-servernum',
                  '--server-num=1',
                  str_c('-s \'-screen 0 ', width, 'x', height, 'x8\''),
                  igv_sh,
                  str_c(igv_args, collapse = ' '),
                  '-b',
                  batch_tmp,
                  sep = ' ')
        # wrap singularity if using
        if (!is.null(singularity_img)) {
            bind_dirs <- c(getwd(), map_chr(c(bams, snapshot_tbl$filename, batch_tmp), dirname)) %>% unique() %>% setdiff('.')
            cmd <-
                str_c(singularity_bin, 
                      'exec',
                      str_c('-B ', bind_dirs, collapse = ' '),
                      singularity_img,
                      cmd,
                      sep = ' ')
        }
        # execute command
        system(cmd)
    }
    
    return(candidates)
}
