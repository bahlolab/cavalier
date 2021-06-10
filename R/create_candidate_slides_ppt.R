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
create_candidate_slides_ppt <- function(candidates, output_dir, output_cols,
                                        genemap2=NULL,
                                        GTEx_median_rpkm=NULL,
                                        GTEx_tissues=NULL,
                                        title_col = NULL,
                                        add_data_col = NULL,
                                        slide_template = NULL) {

    # shouldn't be need as in imports, but getting read_pptx() not found
    if (nrow(candidates) == 0) {
        return(NULL)
    }
    if (is.null(slide_template)) {
        slide_template <- system.file("ppt", "template.pptx", package = "cavalier")
    }
    stopifnot(file.exists(slide_template))
    rownames(candidates) <- 1:nrow(candidates)
    
    if (all(c("reference", "alternate") %in% colnames(candidates))) {
        candidates$"ref / alt" <- paste(candidates$reference, "/", candidates$alternate)
    }
    if ("change" %in% colnames(candidates)) {
        candidates$change <- gsub(";", "\n", candidates$change, fixed=TRUE)
    }
    
    output_cols <- intersect(output_cols, colnames(candidates))

    output_dir <- file.path(output_dir, 'pptx_files')
    dir.create(output_dir, recursive=TRUE, showWarnings = FALSE)
    
    slides <- read_pptx(slide_template)
    
    # caclulate position of slide elements (in inches), based on widescreen slide of 13.33 x 7.5 inches
    tot_width <- 13.33
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
    
    print(slides, target = file.path(output_dir, 'candidate_variants_report.pptx'))
}

