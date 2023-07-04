
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
                                 snapshot_dir = NULL,
                                 vcf_file = NULL,
                                 overwrite = TRUE,
                                 slop = 20,
                                 width = 500,
                                 height = 700,
                                 crop_left = 19,
                                 crop_right = 17,
                                 crop_top = 30,
                                 name_panel_width = 10,
                                 prefs = character()) {
  
  # TODO - check bam files have indexes, as IGV just hangs if they don't
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
    is_scalar_integerish(crop_top),
    is_scalar_integerish(name_panel_width),
    is_character(prefs))
  
  # get cavalier options
  snapshot_dir <- `if`(is.null(snapshot_dir),
                       get_cavalier_opt('snapshot_dir'),
                       snapshot_dir)
  ref_genome <- get_cavalier_opt('ref_genome')
  xvfb_run_cmd <- get_cavalier_opt('xvfb_run_cmd')
  igv_cmd <- get_cavalier_opt('igv_cmd')
  singularity_img <- get_cavalier_opt('singularity_img')
  singularity_cmd <- get_cavalier_opt('singularity_cmd')
  
  # check cavalier options
  assert_that(
    is_scalar_character(snapshot_dir),
    is_scalar_character(ref_genome) && ref_genome %in% c('hg38', 'hg19'),
    is_scalar_character(xvfb_run_cmd),
    is_scalar_character(igv_cmd),
    is_null_or_file(singularity_img),
    is.null(singularity_cmd) | is_scalar_character(singularity_cmd))
  
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
    mutate(filename = file.path(snapshot_dir,
                                str_c(str_c(sample, id, chrom, start, end, sep='_'), '.png')))
  
  to_snap <-
    snapshot_tbl %>% 
    filter(overwrite | !file.exists(filename)) %>% 
    mutate(region = str_c(chrom, ':', start, '-', end)) %>% 
    select(bam, region, filename) %>% 
    distinct() %>% 
    nest(data = -bam) 
  
  if (nrow(to_snap)) {
    
    # create dummy HOME directory to ensure idempotence
    home_dir <- tempfile(tmpdir = snapshot_dir, pattern = 'igv_home_')
    # igv checks for genomes and settings in ~/igv on startup
    igv_dir <- file.path(home_dir, 'igv')
    genome_dir <- file.path(igv_dir, 'genomes')
    assert_create_dir(genome_dir)
    
    # setup IGV genome from cache to avoid downloading repeatedly
    genome_file <- get_igv_genome(ref_genome)
    genome_cp <- file.path(genome_dir, basename(genome_file))
    if (!file.exists(genome_cp)) {
      file.copy(normalizePath(genome_file), genome_dir, overwrite = TRUE)
    }
    
    # width and height before cropping
    tot_width <- width + crop_left + crop_right
    tot_height <- height + crop_top
    
    # write prefs.properties to set IGV dimensions and other defaults
    c(str_c('NAME_PANEL_WIDTH=', name_panel_width),
      str_c('IGV.Bounds=0,0,', tot_width, ',', height),
      str_c('DEFAULT_GENOME_KEY=', ref_genome),
      prefs) %>%
      write_lines(file.path(igv_dir, 'prefs.properties'))
    
    # IGV gets java arguments from ~/.igv/java_arguments
    # here we can set -Duser.home= to change the home directory as detected by IGV
    hidden_igv_dir <- file.path(home_dir, '.igv')
    assert_create_dir(hidden_igv_dir)
    str_c('-Duser.home=', normalizePath(home_dir)) %>% 
      write_lines(file.path(hidden_igv_dir, 'java_arguments'))
    
    batch_file <- file.path(home_dir, 'igv_snapshot.batch')
    
    to_snap %>% 
      pmap(function(bam, data) {
        # IGV batch script
        c('new',
          str_c('load ', c(bam, vcf_file)),
          str_c('maxPanelHeight ', height),
          pmap(data, function(region, filename) {
            c(str_c('goto ', region), str_c('snapshot ', filename))
          }) %>% flatten_chr())
      }) %>% 
      flatten_chr() %>% 
      c(., 'exit') %>% 
      write_lines(batch_file)
    
    # xvfb/ IGV commamnd
    cmd <- 
      str_c(str_c('HOME=',normalizePath(home_dir)),
            ' bash -c "',
            str_c(
              xvfb_run_cmd,
              '--auto-servernum',
              '--server-num=1',
              str_c('-s \'-screen 0 ', tot_width, 'x', height, 'x8\''),
              igv_cmd,
              '-b', batch_file,
              sep = ' '),
            '"')
    
    # wrap singularity is using
    if (!is.null(singularity_img)) {
      # create script to run as we are already wrapping a lot commands
      script <- file.path(home_dir, 'xvfb_run_igv.sh')
      write_lines(cmd, script)
      Sys.chmod(script, '750')
      
      bind_dirs <- 
        c(getwd(),
          snapshot_dir,
          normalizePath(c(vcf_file, bam_files)) %>% dirname()) %>% 
        normalizePath() %>% 
        unique() %>% 
        remove_child_dirs()
      cmd <-
        str_c(singularity_cmd, 
              'exec',
              str_c('-B ', bind_dirs, collapse = ' '),
              # '--home', normalizePath(snapshot_dir),
              singularity_img,
              # cmd,
              script,
              sep = ' ')
    } 
    # execute command
    message('executing: ', cmd)
    assert_that(system(cmd) == 0)
    
    # cleanup dummy home directory
    unlink(home_dir, recursive = TRUE)
  }
  
  # crop images and output
  sample_snaps <-
    snapshot_tbl %>% 
    mutate(cropped = str_replace(filename, '.png', '.cropped.png')) %T>%
    with(crop_pad_png(filename, cropped,
                      crop_left = crop_left,
                      crop_right = crop_right,
                      crop_top = 60,
                      crop_height = height)) %>% 
    select(id, sample, filename = cropped)
}

# download IGV genome file
get_igv_genome <- function(ref_genome) {
  
  cache_dir <- get_cache_dir()
  
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

#' @importFrom ggplot2 ggplot theme_void coord_fixed aes scale_x_continuous scale_y_continuous
#' @importFrom ggimg geom_rect_img
plot_igv_snapshots <- function(sample_snaps,
                               max_cols = 5L,
                               width = 500,
                               height = 700,
                               base_size = 14,
                               strip.position = 'top') {
  # input df of sample, filename
  p <-
    sample_snaps %>% 
    ggplot() + 
    geom_rect_img(aes(xmin = 0, xmax = width, ymin = 0, ymax = height, img = filename)) +
    facet_wrap(~sample, strip.position = strip.position,
               ncol = min(nrow(sample_snaps), max_cols)) + 
    theme_void(base_size = base_size) + 
    coord_fixed() +
    scale_x_continuous(limits = c(0,width), expand = c(0.01, 0.01)) +
    scale_y_continuous(limits = c(0,height), expand = c(0.01, 0.01)) 
  return(p)
}


#' @importFrom png readPNG writePNG
#' @importFrom stringr str_c str_remove
#' @importFrom abind abind
crop_pad_png <- function(input_png,
                         output_png,
                         crop_left = 0, 
                         crop_right = 0,
                         crop_top = 0,
                         crop_height = NULL,
                         pad = TRUE) {
  assert_that(
    is_character(input_png),
    all(file.exists(input_png)),
    is_character(output_png),
    length(input_png) == length(output_png),
    is_scalar_integerish(crop_left),
    is_scalar_integerish(crop_right),
    is_scalar_integerish(crop_top),
    is.null(crop_height) | is_scalar_integerish(crop_right))
  
  walk2(input_png, output_png, function(input, output) {
    png <- readPNG(input)
    # crop left, right and top
    png <- png[seq.int(crop_top +1, dim(png)[1]),
               seq.int(crop_left+1, dim(png)[2] - crop_right),]
    # crop bottom if too tall
    if (!is.null(crop_height) && dim(png)[1] > crop_height) {
      png <- png[1:crop_height, ,]
    }
    # add invisible padding if too short
    if (!is.null(crop_height) && pad && dim(png)[1] < crop_height) {
     png <- 
       abind(png,
             array(1, dim = c(crop_height - dim(png)[1], dim(png)[2], dim(png)[3])),
             along = 1)
    }
    writePNG(png, target = output)
  })
}
