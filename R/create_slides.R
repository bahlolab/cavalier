#' Create powerpoint .pptx output slide for each candidate variant
#' 
#' @param candidates candidate variants data.frame
#' @param output_dir cavalier output directory
#' @param genemap2 location of genemap2.txt file downloaded from OMIM (see https://omim.org/downloads/)
#' @param GTEx_median_rpkm location of GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct.gz file downloaded from GTEx Portal (see https://gtexportal.org/home/datasets)
#' @param GTEx_tissues optionally specify list of tissues to plot GTEx expression data
#' @param hide_missing_igv hide variants that are missing IGV snapshot (default: FALSE)
#' @param layout slide layout choice: "individual" or "multiple" designed for a single or multiple individuals (default: "individual")
#' @importFrom officer read_pptx add_slide ph_with ph_location_type ph_location_template external_img 
#' @importFrom flextable flextable delete_part theme_zebra italic bold colformat_char align autofit fit_to_width
#' @importFrom dplyr arrange select mutate
#' @importFrom tibble rownames_to_column
#' @importFrom purrr walk map
create_slides <- function(title_col,
                          output = 'cavalier_slides.pptx',
                          bam_files = NULL,
                          ped_file = NULL,
                          vcf_file = NULL,
                          genemap2_file = NULL,
                          layout = layout_single(),
                          title_col = 'title',
                          slide_template = slide_template(),
                          var_info = var_info()) 
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
    
    elements <- map(layout$slides, ~ map(.$rows, 'elements')) %>% unlist()
    
    # check we have required data for plot elements
    assert_that(
        ! 'igv' %in% elements | !is.null(bam_files),
        ! 'omim' %in% elements | !is.null(genemap2_file),
        ! 'pedigree' %in% elements | !is.null(ped_file))
    
    slide_data <- tibble(id = seq_len(nrow(variants)),
                         title = variants[[title_col]])
    
    if ('var_info' %in% elements) {
        
        names(var_info) <- 
            `if`(is.null(names(var_info)),
                 var_info,
                 if_else(names(var_info) == '', var_info, names(var_info)))
        
        slide_data <-
            select(variants, all_of(var_info)) %>% 
            mutate(id = slide_data$id) %>% 
            nest(var_info = -id) %>% 
            full_join(slide_data, by = 'id')
    }
    
    if ('igv' %in% elements) {
        slide_data <-
            create_igv_snapshots(variants, bam_files) %>% 
            mutate(id = slide_data$id) %>% 
            nest(igv = -id) %>% 
            full_join(slide_data, by = 'id')
    }
    
    if ('pedigree' %in% elements) {
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
                    plot_ped(draw = FALSE)
            })
    }
    
    if ('gtex' %in% elements) {
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
    
    if ('omim' %in% elements) {
        slide_data <-
            variants %>% 
            select(gene) %>%
            mutate(id = seq_along(gene)) %>% 
            left_join(get_omim_genemap2(genemap2_file) %>% 
                          select(gene = symbol, phenotype, inheritance),
                      by = 'gene') %>% 
            select(-gene) %>% 
            nest(omim = -id) %>% 
            full_join(slide_data, by = 'id')
    }
    
    custom <- elements[str_starts(elements, 'custom_')]
    if (length(custom)) {
        sel <- setNames(str_remove(custom, '^custom_'), custom)
        custom_data <-
            select(variants, all_of(sel)) %T>%
            # check correct data types
            (function(data) {
                assert_that(
                    all(map_lgl(data, is.list)),
                    all(map_lgl(data, ~ all(map_lgl(., ~ {
                        is.data.frame(.) || is(., 'grob') || is (., 'gg')
                    }))))
                )
            })
        slide_data <- bind_cols(slide_data, custom_data)
    }
    
    slides <- read_pptx(slide_template)
    
    slide_data %>% 
        pwalk(function(...) {
            add_slides(slides, layout, dots_list(...))
        }) 
    
    print(slides, target = file.path(output_dir, 'candidate_variants_report.pptx'))
    
    return(invisible(variants))
}

# add slides using layout
add_slides <- function(slides, layout, data) 
{
    # caclulate position of slide elements (in inches), based on widescreen slide of 13.33 x 7.5 inches
    tot_width <- 13.333
    tot_height <- 7.5
    pad <- 0.02
    # heights
    pad_h <- tot_height * pad
    header <- 0.8
    r1_p <- 2/3
    nrows <- 2
    use_height <- tot_height - header - (nrows+1) * pad_h
    r1_h <- r1_p * use_height
    r2_h <- use_height - r1_h
    # widths, top row
    pad_w <- pad * tot_width
    n_cols_1 <- 3
    use_width_1 <- tot_width - (n_cols_1+1) * pad_w
    seg_width_1 <- use_width_1 / n_cols_1
    seg_left_1 <- seq(from=pad_w, by = seg_width_1 + pad_w, length.out = n_cols_1)
    # widths, bottom row
    n_cols_2 <- 2
    omim_p <- 2/5
    use_width_2 <- tot_width - (n_cols_2+1)* pad_w
    omim_width <- use_width_2 * omim_p
    add_data_width <- use_width_2 - omim_width
    
    candidates %>% 
        arrange(`inheritance model`, gene) %>% 
        split.data.frame(seq_len(nrow(.))) %>% 
        walk(function(cand) {
            gene <- cand$gene
            # Information table
            info_ft <-
                select(cand, all_of(output_cols)) %>% 
                t() %>% 
                as.data.frame(stringsAsFactors = FALSE) %>% 
                rownames_to_column('name') %>% 
                mutate_if(is.character, ~str_replace_all(., '\\n', ', ')) %>% 
                flextable() %>% 
                delete_part(part = "header") %>% 
                theme_zebra(even_header = 'white', even_body = 'white') %>% 
                italic(j = 1) %>% 
                bold(j = 1) %>% 
                colformat_char(j = 1, suffix = ':') %>% 
                align(j = 1, align = 'right', part = 'all') %>% 
                flextable_fit(width = seg_width_1, height = r1_h)
            
            omim_ft <- NULL
            omim_table <- omim_table(gene, genemap2=genemap2, wrap = FALSE)
            if (!is.null(omim_table)) {
                omim_ft <-
                    omim_table %>% 
                    flextable() %>% 
                    theme_zebra(even_header = 'white', even_body = 'white') %>% 
                    flextable_fit(width = omim_width, 
                                  height = r2_h)
            }
            
            add_data_ft <- NULL
            if (!is.null(add_data_col)) {
                add_data_table <- cand[[add_data_col]][[1]]
                if (nrow(add_data_table)) {
                    add_data_ft <-
                        add_data_table %>% 
                        flextable() %>% 
                        theme_zebra(even_header = 'white', even_body = 'white') %>% 
                        autofit() %>% 
                        flextable_fit(width = add_data_width, 
                                      height = r2_h)
                }
            }
            
            gtex_plot <- plot_gtex_expression(gene, GTEx_median_rpkm=GTEx_median_rpkm, tissues=GTEx_tissues)
            
            if (!is.null(title_col)) {
                title <- cand[[title_col]][[1]]
            } else {
                title <- gene
            }
            
            # calculate igv height to preserve aspect ratio
            igv_img <- read_png(cand$igv_filename)
            igv_height <- with(attributes(igv_img)$dims, height * (seg_width_1 / width))
            
            slides <-
                slides %>% 
                add_slide(layout = "Title and Content") %>% 
                ph_with(value = title,
                        location = ph_location_type(type = "title")) %>% 
                ph_with(value = info_ft,
                        location = ph_location_template(left = seg_left_1[1],
                                                        top = header + pad_h,
                                                        width = seg_width_1,
                                                        height = r1_h)) %>% 
                ph_with(value = igv_img,
                        location = ph_location_template(left = seg_left_1[2],
                                                        top = header + pad_h,
                                                        width = seg_width_1,
                                                        height = igv_height)) %>% 
                { `if`(is.null(gtex_plot), .,
                       ph_with(.,
                               value = gtex_plot,
                               location = ph_location_template(left = seg_left_1[3],
                                                               top = header + pad_h,
                                                               width = seg_width_1,
                                                               height = r1_h))
                )} %>% 
                { `if`(is.null(omim_ft), .,
                       ph_with(.,
                               value = omim_ft,
                               location = ph_location_template(left = seg_left_1[1],
                                                               top = header + r1_h + 2*pad_h,
                                                               width = omim_width))
                )} %>% 
                { `if`(is.null(add_data_ft), .,
                       ph_with(.,
                               value = add_data_ft,
                               location = ph_location_template(left = seg_left_1[1] + omim_width + pad_w,
                                                               top = header + r1_h + 2*pad_h,
                                                               width = add_data_width))
                )} 
        })
    
}

#' @export
slide_template <- function() {
    system.file("ppt", "template.pptx", package = "cavalier")
}

#' @export
layout_single <- function(extra = NULL,
                          ...)
{
    assert_that(is.null(extra) | (is_scalar_character(extra) && is_valid_slide_element(extra)))
    
    slide_layout(
        slide(c('var_info', 'igv', 'gtex'),
              c('omim', extra),
              heights = c(2,1)),
        ...
    )
}

#' @export
layout_multiple <- function(extra = NULL,
                            pedigree = FALSE,
                            ...) 
{
    assert_that(is.null(extra) | (is_scalar_character(extra) && is_valid_slide_element(extra)))
    
    slide_layout(
        slide(c('var_info', `if`(pedigree, 'pedigree', NULL), 'gtex'),
              c('omim', extra),
              heights = c(2,1)),
        slide('igv'),
        ...
    )
}

#' @importFrom rlang is_scalar_double dots_list
#' @importFrom purrr map_lgl
#' @export
slide_layout <- function(...,
                         slide_height = 7.5,
                         slide_width = 13.333,
                         pad = 0.02)
{
    slides <- dots_list(...)
    assert_that(all(map_lgl(slides, ~ is(., 'slide'))),
                is_scalar_double(slide_height),
                is_scalar_double(slide_width),
                is_scalar_double(pad),
                pad >= 0,
                pad < 1)
    layout <- as.list(environment())
    class(layout) <- c('list', 'layout')
    return(layout)
}

#' @importFrom rlang dots_list
#' @importFrom purrr map_lgl walk
#' @importFrom stringr str_starts
#' @export
slide <- function(...,
                  heights = NULL)
{
    rows <- dots_list(...)
    
    assert_that(length(rows) > 0,
                all(map_lgl(rows, is_valid_row)),
                is.null(heights) | (is_number(heights) & length(heights) == length(rows)))
    
    heights <- `if`(is.null(heights),
                    rep(1/length(rows), length(rows)),
                    heights / sum(heights))
    
    rows <- map(rows, function(row) {
        `if`(is_named(row), 
             list(elements = names(row),
                  widths = row / sum(row)),
             list(elements = row,
                  widths = rep(1/length(row), length(row))))
    })
    
    walk(rows, function(row) {
        assert_that(all(is_valid_slide_element(row$elements)))
    })
    
    sl <- list(rows = rows, heights = heights)
    class(sl) <- c('list', 'slide')
    return(sl)
}

slide_elements <- function() {
    c('igv', 'var_info', 'gtex', 'omim', 'pedigree') 
}

is_valid_slide_element <- function(x) {
    x %in% slide_elements() | str_starts(x, 'custom_')
}

#' @importFrom rlang is_character is_double is_integer is_named
is_valid_row <- function(x) {
    is_character(x) | ( is_number(x) & is_named(x) )
}


