
#' @importFrom officer read_pptx add_slide ph_with ph_location_type ph_location_template external_img 
#' @importFrom dplyr arrange select mutate
#' @importFrom tibble rownames_to_column
#' @importFrom purrr walk map
#' @export
create_slides <- function(variants,
                          output = 'cavalier_slides.pptx',
                          bam_files = NULL,
                          ped_file = NULL,
                          vcf_file = NULL,
                          layout = layout_single(),
                          title_col = 'title',
                          slide_template = get_slide_template(),
                          var_info = get_var_info())
{
  # check args
  assert_that(
    is.data.frame(variants),
    is_null_or_files(bam_files, named = TRUE),
    is_null_or_file(ped_file),
    is_null_or_file(vcf_file),
    is_scalar_character(output),
    is_scalar_character(slide_template) && file.exists(slide_template),
    is_scalar_character(title_col),
    is_character(var_info))
  
  # check we have required data for plot elements
  assert_that(
    ! 'igv' %in% layout$element | !is.null(bam_files),
    ! 'omim' %in% layout$element | !is.null(genemap2_file),
    ! 'pedigree' %in% layout$element | !is.null(ped_file))
  
  slide_data <- tibble(id = seq_len(nrow(variants)),
                       title = variants[[title_col]])
  
  if ('var_info' %in% layout$element) {
    
    names(var_info) <- 
      `if`(is.null(names(var_info)),
           var_info,
           if_else(names(var_info) == '', 
                   var_info, names(var_info)))
    
    slide_data$var_info <-
      select(variants, all_of(var_info)) %>% 
      (function(x) {
        # add in gnomad link
        `if`('af_gnomad' %in% var_info,
             names(var_info)[which(var_info == 'af_gnomad')] %>% 
               str_c( '_url') %>% 
               { mutate(x, !!. := gnomad_link(variants)) },
             x)
      }) %>% 
      (function(x) {
        # add in dbsnp link
        `if`('db_snp' %in% var_info,
             names(var_info)[which(var_info == 'db_snp')] %>% 
               str_c( '_url') %>% 
               { mutate(x, !!. := dbsnp_link(variants$db_snp)) },
             x)
      }) %>% 
      (function(x) {
        # add in genecards
        `if`('gene' %in% var_info,
             names(var_info)[which(var_info == 'gene')] %>% 
               str_c( '_url') %>% 
               { mutate(x, !!. := genecards_link(variants$gene)) },
             x)
      }) %>% 
      (function(x) {
        # add in genecards
        `if`('ensembl_gene' %in% var_info,
             names(var_info)[which(var_info == 'ensembl_gene')] %>% 
               str_c( '_url') %>% 
               { mutate(x, !!. := ensembl_gene_link(variants$ensembl_gene)) },
             x)
      }) %>% 
      mutate(id = slide_data$id) %>% 
      nest(var_info = -id) %>% 
      with(map(var_info, flex_table, transpose = TRUE))
  }
  
  if ('igv' %in% layout$element) {
    
    slide_data$igv <-
      create_igv_snapshots(variants, bam_files,
                           width = 500,
                           height = 700) %>% 
      nest(data = -id) %>% 
      mutate(data = map2(id, data, function(id, data) {
        variants$genotype[id, ] %>% 
          pivot_longer(everything(),
                       names_to = 'sample',
                       values_to = 'gt') %>% 
          right_join(data, by = 'sample') %>% 
          mutate(sample = str_c(sample, ': ', gt))
      })) %>% 
      with(map(data, plot_igv_snapshots, width = 500, height = 700))
  }
  
  if ('pedigree' %in% layout$element) {
    
    ped_df <- read_ped(ped_file)
    slide_data$pedigree <-
      seq_len(nrow(variants)) %>% 
      map(function(i) {
        variants$genotype[i, ] %>% 
          pivot_longer(everything(),
                       names_to = 'iid',
                       values_to = 'gt') %>% 
          right_join(ped_df, by = 'iid') %>% 
          mutate(gt = replace_na(gt, 'ND'),
                 label = str_c(iid, gt, sep = '\n')) %>% 
          plot_ped()
      })
  }
  
  if ('gtex' %in% layout$element) {
    
    slide_data$gtex <-
      variants %>% 
      select(ensembl_gene, gene) %>% 
      pmap(function(ensembl_gene, gene) {
        if (!is.na(ensembl_gene)) {
          plot_gtex_expression(gene, ensembl_id = ensembl_gene)
        } else if(!is.na(gene)) {
          plot_gtex_expression(gene,)
        } else {
          ggdraw()
        }
      })
  }
  
  if ('omim' %in% layout$element) {
    
    slide_data$omim <-
      variants %>%
      select(gene) %>%
      mutate(id = slide_data$id) %>%
      left_join(get_omim_gene_map(), by = c(gene = 'symbol')) %>% 
      mutate(disease_name = disease_names(disease_id)) %>% 
      mutate(disease_name_url = str_c('https://omim.org/entry/', str_extract(disease_id, '\\d+$')),
             disease_id_url = disease_name_url) %>% 
      select(id, disease_id, disease_name, inheritance, disease_id_url, disease_name_url) %>% 
      nest(omim = -id) %>%
      with(map(omim, flex_table))
  }
  
  custom <- layout$element[str_starts(layout$element, 'custom_')]
  if (length(custom)) {
    
    sel <- setNames(str_remove(custom, '^custom_'), custom)
    custom_data <-
      select(variants, all_of(sel))
    # %T>%
    # # check correct data types
    # (function(data) {
    #     assert_that(
    #         all(map_lgl(data, is.list)),
    #         all(map_lgl(data, ~ all(map_lgl(., ~ {
    #             is.data.frame(.) || is (., 'gg') | is(., 'flextable')
    #         }))))
    #     )
    # })
    slide_data <- bind_cols(slide_data, custom_data)
  }
  
  slides <- read_pptx(slide_template)
  
  slide_data %>% 
    pwalk(function(...) {
      # data <<- dots_list(...)
      slides <- add_slides(slides, layout, dots_list(...))
    })
  
  if (nrow(slide_data) == 0) {
    # add slide stating no results
    slides <-
      slides %>% 
      add_slide(layout = "Title and Content") %>% 
      ph_with(value = 'No Variants Found',
              location = ph_location_type(type = "title"))
  }
  
  print(slides, target = output)
  re_encode_pptx_hlinks(output)
  
  return(invisible(variants))
}

# add slides using layout
add_slides <- function(slides, layout, data)
{
  n_slides <- n_distinct(layout$slide_num)
  
  title <- 
    `if`(n_slides  == 1,
         data$title,
         str_c(data$title, ' (', seq_len(n_slides), '/', n_slides, ')'))
  
  walk(seq_len(n_slides),  function(i) {
    # add title
    slides <-
      slides %>% 
      add_slide(layout = "Title and Content") %>% 
      ph_with(value = title[i],
              location = ph_location_type(type = "title"))
    
    layout %>%
      filter(slide_num == i) %>%
      pwalk(function(element, x_left, width, y_top, height, transpose, ...) {
        value <- data[[element]]
        
        if (is.data.frame(value)) {
          value <- flex_table(value)
        }
        if (is(value, 'flextable')) {
          value <- fit_flex_table(value, width = width, height = height,
                                  expand_rows = element == 'var_info')
        }
        # add item to slides
        if (!is.null(value)) {
          slides <-
            slides %>% 
            ph_with(value = value,
                    location = ph_location_template(
                      left = x_left,
                      top = y_top,
                      width = width,
                      height = height))
        }
      })
  })
  
  return(slides)
}

#' @export
get_slide_template <- function() {
  system.file("ppt", "template.pptx", package = "cavalier")
}

#' @export
layout_single <- function(...)
{
  slide_layout(
    c('var_info', 'igv', 'gtex'),
    c('omim'),
    heights = c(2,1),
    ...)
}

#' @export
layout_multiple <- function(pedigree = FALSE,
                            ...) 
{
  bind_rows(
    slide_layout(c('var_info', `if`(pedigree, 'pedigree', NULL), 'gtex'),
                 c('omim'),
                 heights = c(2,1),
                 ...),
    slide_layout(c('igv'),
                 slide_num = 2L,
                 ...))
}

#' @importFrom rlang dots_list is_scalar_double
#' @importFrom purrr map_lgl walk map_df
#' @importFrom stringr str_starts
#' @export
slide_layout <- function(...,
                         heights = NULL,
                         title_height = 0.1,
                         slide_height = 7.5,
                         slide_width = 13.333,
                         pad = 0.02,
                         slide_num = 1L,
                         transpose = 'var_info')
{
  rows <- dots_list(...)
  
  assert_that(
    is_scalar_double(title_height), title_height >= 0, title_height < 1,
    is_scalar_double(slide_height),
    is_scalar_double(slide_width),
    is_scalar_double(pad), pad >= 0, pad < 1,
    length(rows) > 0,
    all(map_lgl(rows, is_valid_row)),
    is.null(heights) | (is_number(heights) & length(heights) == length(rows)))
  
  
  # get coordinates for each element
  layout_df <-
    map_df(rows, function(row) {
      `if`(is_named(row), 
           tibble(element = names(row),
                  width = row / sum(row)),
           tibble(element = row,
                  width = rep(1/length(row), length(row)))) %>% 
        pad_x(pad = pad * slide_width, 
              slide_width = slide_width) %>% 
        nest(data = everything())
    }) %>% 
    mutate(height = `if`(is.null(heights),
                         rep(1/length(rows), length(rows)),
                         heights / sum(heights))) %>% 
    pad_y(pad = pad * slide_width, 
          title_height = title_height * slide_height,
          slide_height = slide_height) %>% 
    unnest(data) %>% 
    mutate(transpose = element %in% transpose,
           slide_num = slide_num)
  
  return(layout_df)
}

pad_x <- function(col_df, pad, slide_width) {
  mutate(col_df, 
         abs = FALSE,
         row_num = row_number() * 2L) %>% 
    bind_rows(tibble(element = 'pad',
                     width = pad,
                     abs = TRUE,
                     row_num = seq.int(1, nrow(col_df) * 2 +1, by = 2))) %>% 
    arrange(row_num) %>% 
    mutate(width = if_else(abs, width, width * (slide_width - sum(width[abs]))),
           x_right = cumsum(width),
           x_left = x_right - width,
           x_cen = (x_right + x_left) / 2) %>% 
    select(element, width, x_left, x_cen, x_right) %>% 
    filter(element != 'pad')
}

pad_y <- function(row_df, pad, title_height, slide_height) {
  mutate(row_df, 
         abs = FALSE,
         row_num = row_number() * 2L) %>% 
    bind_rows(tibble(data = list(NULL),
                     height = pad,
                     abs = TRUE,
                     row_num = seq.int(1, nrow(row_df) * 2 +1, by = 2))) %>% 
    add_row(row_num = 0L,
            data = NULL, 
            abs = TRUE,
            height = title_height) %>% 
    arrange(row_num) %>% 
    mutate(height = if_else(abs, height, height * (slide_height - sum(height[abs]))),
           y_bot = cumsum(height),
           y_top = y_bot - height,
           y_cen = (y_bot + y_top) / 2) %>%
    select(data, height, y_bot, y_cen, y_top) %>%
    filter(map_lgl(data, ~ ! is.null(.)))
}



slide_elements <- function() {
  c('igv', 'var_info', 'gtex', 'omim', 'pedigree') 
}

is_valid_slide_element <- function(x) {
  x %in% slide_elements() | str_starts(x, 'custom_')
}

#' @importFrom rlang is_character is_double is_integer is_named
is_valid_row <- function(x) {
  (is_character(x) & all(is_valid_slide_element(x))) | 
    ( is_number(x) && is_named(x) && all(is_valid_slide_element(names(x))))
}

#' @export
get_var_info <- function() {
  c(`Gene Symbol` = 'gene',
    `Ensembl Gene` = 'ensembl_gene',
    Inheritance = 'inheritance',
    Consequence = 'consequence',
    dbSNP = 'db_snp',
    HGVSg = 'hgvs_genomic', 
    HGVSc = 'hgvs_coding',
    HGVSp = 'hgvs_protein',
    SIFT = 'sift',
    PolyPhen = 'polyphen',
    `gnomAD AF` = 'af_gnomad'
  )    
}

re_encode_pptx_hlinks <- function(target) 
{
  hlink_type <- "http://schemas.openxmlformats.org/officeDocument/2006/relationships/hyperlink"
  tmp_dir <- tempfile(pattern = '.', tmpdir = '.')
  zip::unzip(target, exdir = tmp_dir)
  
  changed <- FALSE
  
  list.files(file.path(tmp_dir, 'ppt/slides/_rels'),
             full.names = TRUE) %>% 
    walk(function(file) {
      content <- 
        xml2::read_xml(file) %>% 
        xml2::as_list()
      
      rep_rel <-
        map(content$Relationships, function(item) {
          atrib <- attributes(item)
          if (atrib$Type == hlink_type) {
            if (str_detect(atrib$Target, '%')) {
              atrib$Target <- 
                utils::URLdecode(atrib$Target) %>% 
                utils::URLencode()
              changed <<- TRUE
            }
          }
          attributes(item) <- atrib
          item
        })
      
      attributes(rep_rel) <- attributes(content$Relationships)
      content$Relationships <- rep_rel
      
      if (changed) {
        xml2::write_xml(xml2::as_xml_document(content), file)
      }
    })
  
  if (changed) {
    officer::pack_folder(tmp_dir, target)
  }
  unlink(tmp_dir, recursive = TRUE)
  
  invisible(NULL)
}

