
# creates IGV snapshots using xvfb-run and igv.sh
# Can either set igv_sh to full path to igv.sh, or Sys.setenv(PATH=...), likewise for xvfb-run/ singularity
#' @importFrom purrr walk pwalk pmap flatten_chr map_chr 
#' @importFrom stringr str_c
#' @importFrom dplyr distinct transmute mutate filter
#' @importFrom tidyr expand_grid pivot_wider
#' @importFrom readr write_lines
#' @importFrom magrittr '%T>%'
#' @export
create_igv_snapshots <- function(variants, bam_files,
                                 vcf_file = NULL,
                                 overwrite = FALSE,
                                 slop = 20,
                                 width = 720,
                                 height = 620,
                                 crop_left = 20,
                                 crop_right = 17,
                                 name_panel_width = 10,
                                 prefs = character()) {
    
    # TODO - store IGV genome files in cache_dir, pre-download from broad link before running IGV
    # check args
    assert_that(
        is.data.frame(variants),
        is_null_or_files(bam_files, named = TRUE),
        is_null_or_file(vcf_file),
        is_bool(overwrite),
        is_scalar_integerish(slop),
        is_scalar_integerish(width),
        is_scalar_integerish(height),
        is_scalar_integerish(crop_left),
        is_scalar_integerish(crop_right),
        is_scalar_integerish(name_panel_width),
        is_character(prefs))
        
    # get cavalier options
    igv_snapshot_dir <- get_cavalier_opt('igv_snapshot_dir')
    ref_genome <- get_cavalier_opt('ref_genome')
    xvfb_run_cmd <- get_cavalier_opt('xvfb_run_cmd')
    igv_cmd <- get_cavalier_opt('igv_cmd')
    singularity_img <- get_cavalier_opt('singularity_img')
    singularity_cmd <- get_cavalier_opt('singularity_cmd')
    
    # check cavalier options
    assert_that(
        is_scalar_character(igv_snapshot_dir),
        is_scalar_character(ref_genome) && ref_genome %in% c('hg38', 'hg19'),
        is_scalar_character(xvfb_run_cmd),
        is_scalar_character(igv_cmd),
        is_null_or_file(singularity_img),
        is.null(singularity_cmd) | is_scalar_character(singularity_cmd))
    
    # setup working directories
    user_dir <- file.path(igv_snapshot_dir, 'igv')
    genome_dir <- file.path(user_dir, 'genomes')
    assert_create_dir(genome_dir, recursive = TRUE)
    # setup IGV genome from cache to avoid downloading 
    genome_file <- get_igv_genome(ref_genome)
    genome_link <- file.path(genome_dir, basename(genome_file))
    if (!file.exists(genome_link)) {
      file.copy(normalizePath(genome_file), genome_dir, overwrite = TRUE)
    }
    
    # write prefs.properties to set IGV dimensions
    c(str_c('NAME_PANEL_WIDTH=', name_panel_width),
      str_c('IGV.Bounds=0,0,', width, ',', height),
      str_c('DEFAULT_GENOME_KEY=', ref_genome),
      prefs) %>%
        write_lines(file.path(user_dir, 'prefs.properties'))
    
    snapshot_tbl <-
        variants %>% 
        transmute(id = seq_along(chrom),
                  chrom = chrom,
                  start = pos,
                  end = start + nchar(ref)) %>% 
        mutate(start = start - slop,
               end = end + slop) %>% 
        expand_grid(tibble(sample = names(bam_files),
                           bam = bam_files)) %>% 
        mutate(filename = file.path(igv_snapshot_dir,
                                    str_c(str_c(sample, chrom, start, end, sep='_'), '.png')))
    
    to_snap <-
      snapshot_tbl %>% 
      filter(overwrite | !file.exists(filename)) %>% 
      mutate(region = str_c(chrom, ':', start, '-', end)) %>% 
      select(bam, region, filename) %>% 
      distinct() %>% 
      nest(data = -bam) 
    
    if (nrow(to_snap)) {
      
      batch_file <- file.path(igv_snapshot_dir, 'igv_snapshot.batch')
      
      to_snap %>% 
        pmap(function(bam, data) {
          # IGV batch script
          c('new',
            str_c('load ', c(bam, vcf_file)),
            str_c('maxPanelHeight ', height - 95),
            pmap(data, function(region, filename) {
              c(str_c('goto ', region), str_c('snapshot ', filename))
            }) %>% flatten_chr())
        }) %>% 
        flatten_chr() %>% 
        c(., 'exit') %>% 
        write_lines(batch_file)
      
      # xvfb/ IGV commamnd
      cmd <- 
        str_c(xvfb_run_cmd,
              '--auto-servernum',
              '--server-num=1',
              str_c('-s \'-screen 0 ', width, 'x', height, 'x8\''),
              igv_cmd,
              '--igvDirectory', user_dir,
              '-b', batch_file,
              sep = ' ')
      
      # wrap singularity is using
      if (!is.null(singularity_img)) {
        bind_dirs <- 
          c(getwd(),
            igv_snapshot_dir,
            normalizePath(c(vcf_file, bam_files)) %>% dirname()) %>% 
          normalizePath() %>% 
          unique() %>% 
          remove_child_dirs()
        cmd <-
          str_c(singularity_cmd, 
                'exec',
                str_c('-B ', bind_dirs, collapse = ' '),
                singularity_img,
                cmd,
                sep = ' ')
      } 
      # execute command
      message('executing: ', cmd)
      assert_that(system(cmd) == 0)
    }
    
    # crop images
    if (crop_left > 0 | crop_right > 0) {
      snapshot_tbl <-
        snapshot_tbl %>% 
        mutate(cropped = str_replace(filename, '.png', '.cropped.png')) %T>%
        with(crop_png(filename, cropped)) %>% 
        select(id, sample, filename = cropped)
    }
    
    snapshot_tbl %>% 
      select(id, sample, filename)
}

# download IGV genome file
get_igv_genome <- function(ref_genome) {
  
  cache_dir <- get_cavalier_opt('cache_dir')
  
  genome_uri <- `if`(ref_genome == 'hg38',
                     get_cavalier_opt('igv_hg38_uri'),
                     get_cavalier_opt('igv_hg19_uri'))
  
  genome_file <- file.path(cache_dir, basename(genome_uri))
  
  if (!file.exists(genome_file)) {
    tmp <- tempfile(tmpdir = cache_dir)
    download.file(url = genome_uri, destfile = tmp)
    file.rename(tmp, genome_file)
  }
  
  assert_that(file.exists(genome_file))
  genome_file
}
#' @importFrom png readPNG 
#' @importFrom grid rasterGrob 
#' @importFrom cowplot plot_grid 
arrange_igv_snapshots <- function(sample_pngs,
                                  label = TRUE,
                                  max_cols = 5L) {
  # input df of sample, filename
  p <-
    sample_pngs %>% 
    mutate(plots = map2(sample, filename, function(sm, fn) {
      g <- rasterGrob(readPNG(fn), interpolate=TRUE, just = 'top')
      `if`(label,
           plot_grid(ggdraw() + draw_text(sm), g,
                     ncol = 1, rel_heights = c(1,9)),
           g)
    })) %>% 
    with(plot_grid(plotlist = plots,
                   ncol = min(max_cols, length(plots))))
  return(p)
}


#' @importFrom png readPNG writePNG
#' @importFrom stringr str_c str_remove
crop_png <- function(input_png,
                     output_png,
                     left = 0, 
                     right = 0) {
  assert_that(
    is_character(input_png),
    all(file.exists(input_png)),
    is_character(output_png),
    length(input_png) == length(output_png),
    is_scalar_integerish(left),
    is_scalar_integerish(right))
  
  walk2(input_png, output_png, function(input, output) {
    png <- readPNG(input)
    writePNG(png[,seq.int(left+1, dim(png)[2] - right),],
            target = output)
  })
}
