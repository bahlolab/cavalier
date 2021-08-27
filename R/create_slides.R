
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
                          genemap2_file = NULL,
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
        is_null_or_file(genemap2_file),
        is_scalar_character(output),
        is_scalar_character(slide_template) && file.exists(slide_template),
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
                 if_else(names(var_info) == '', var_info, names(var_info)))
        
        slide_data$var_info <-
            select(variants, all_of(var_info)) %>% 
            mutate(id = slide_data$id) %>% 
            nest(var_info = -id) %>% 
            with(map(var_info, flextable_trans))
    }
    
    if ('igv' %in% layout$element) {
        
        slide_data$igv <-
            create_igv_snapshots(variants, bam_files,
                                 width = 500,
                                 height = 800) %>% 
            nest(data = -id) %>% 
            mutate(data = map2(id, data, function(id, data) {
                variants$genotype[id, ] %>% 
                    pivot_longer(everything(),
                                 names_to = 'sample',
                                 values_to = 'gt') %>% 
                    right_join(data, by = 'sample') %>% 
                    mutate(sample = str_c(sample, ': ', gt))
            })) %>% 
            with(map(data, arrange_igv_snapshots))
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
            left_join(get_omim_genemap2(genemap2_file) %>% 
                          select(gene = symbol,
                                 # mim_number,
                                 OMIM_Phenotype = phenotype,
                                 OMIM_Inheritance = inheritance),
                      by = 'gene') %>% 
            select(-gene) %>% 
            # mutate(url = str_c('https://omim.org/entry/', mim_number)) %>% 
            nest(omim = -id) %>% 
            pull(omim)
            # with(map(omim, function(omim) {
            #     omim %>% 
            #         flextable_reg(col_keys = c('OMIM_Phenotype',
            #                                    'OMIM_Inheritance')) %>% 
            #         compose(j = 'OMIM_Phenotype',
            #                 value = as_paragraph(
            #                     hyperlink_text(x = OMIM_Phenotype, 
            #                                    url = url,
            #                                    props = officer::fp_text(
            #                                        color = 'blue',
            #                                        underlined = TRUE
            #                                    ))))
            # }))
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
                    value <- flextable_reg(value)
                }
                if (is(value, 'flextable')) {
                    value <- flextable_fit(value, width = width, height = height)
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
    c(Gene = 'gene', Inheritance = 'inheritance', Consequence = 'consequence',
      dbSNP = 'db_snp', HGVSg = 'hgvs_genomic', HGVSc = 'hgvs_coding', HGVSp = 'hgvs_protein',
      Grantham = 'grantham_score', SIFT = 'sift', PolyPhen = 'polyphen', RVIS = 'rvis_percentile',
      gnomAD_AF = 'af_gnomad'
    )    
}



