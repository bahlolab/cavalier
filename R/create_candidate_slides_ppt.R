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
    
    # caclulate position of slide elements
    tot_width <- 13.33
    pad <- 0.025
    pad_w <- pad * tot_width
    n_seg <- 3
    use_width <- tot_width - (n_seg+1) * pad_w
    seg_width <- use_width / n_seg
    seg_left <- seq(from=pad_w, by = seg_width + pad_w, length.out = n_seg)
    h1 <- 4
    top <- 0.9
    pad_h <- pad_w
    add_data_width <- 7.2
    omim_width <- 5
    h2 <- 2.1
    
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
                theme_zebra() %>% 
                italic(j = 1) %>% 
                bold(j = 1) %>% 
                colformat_char(j = 1, suffix = ':') %>% 
                align(j = 1, align = 'right', part = 'all') %>% 
                flextable_fit(width = seg_width, height = h1)
            
            omim_ft <- NULL
            omim_table <- omim_table(gene, genemap2=genemap2, wrap = FALSE)
            if (!is.null(omim_table)) {
                omim_ft <-
                    omim_table %>% 
                    flextable() %>% 
                    theme_zebra() %>% 
                    flextable_fit(width = omim_width, 
                                  height = h2)
            }
            
            add_data_ft <- NULL
            if (!is.null(add_data_col)) {
                add_data_table <- cand[[add_data_col]][[1]]
                if (nrow(add_data_table)) {
                    add_data_ft <-
                        add_data_table %>% 
                        flextable() %>% 
                        theme_zebra() %>% 
                        autofit() %>% 
                        flextable_fit(width = add_data_width, 
                                      height = h2)
                }
            }
            
            gtex_plot <- plot_gtex_expression(gene, GTEx_median_rpkm=GTEx_median_rpkm, tissues=GTEx_tissues)
            
            if (!is.null(title_col)) {
                title <- cand[[title_col]][[1]]
            } else {
                title <- gene
            }
            
            slides <-
                slides %>% 
                add_slide(layout = "Title and Content") %>% 
                ph_with(value = title,
                        location = ph_location_type(type = "title")) %>% 
                ph_with(value = info_ft,
                        location = ph_location_template(left = seg_left[1],
                                                        top = top,
                                                        width = seg_width,
                                                        height = h1)) %>% 
                ph_with(value = external_img(cand$igv_filename),
                        location = ph_location_template(left = seg_left[2],
                                                        top = top,
                                                        width = seg_width*1.1)) %>% 
                { `if`(is.null(gtex_plot), .,
                       ph_with(.,
                               value = gtex_plot,
                               location = ph_location_template(left = seg_left[3],
                                                               top = top,
                                                               width = seg_width,
                                                               height = h1))
                )} %>% 
                { `if`(is.null(omim_ft), .,
                       ph_with(.,
                               value = omim_ft,
                               location = ph_location_template(left = seg_left[1],
                                                               top = top + h1 + pad_h,
                                                               width = omim_width))
                )} %>% 
                { `if`(is.null(add_data_ft), .,
                       ph_with(.,
                               value = add_data_ft,
                               location = ph_location_template(left = seg_left[1] + omim_width + 0.35,
                                                               top = top + h1 + pad_h,
                                                               width = add_data_width))
                )} 
        })
    
    slides %>% print(target = file.path(output_dir, 'candidate_variants_report.pptx'))
}

