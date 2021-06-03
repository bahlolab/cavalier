#' Create PDF output slide for each candidate variant
#' 
#' @param candidates candidate variants data.frame
#' @param output_dir cavalier output directory
#' @param genemap2 location of genemap2.txt file downloaded from OMIM (see https://omim.org/downloads/)
#' @param GTEx_median_rpkm location of GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct.gz file downloaded from GTEx Portal (see https://gtexportal.org/home/datasets)
#' @param GTEx_tissues optionally specify list of tissues to plot GTEx expression data
#' @param hide_missing_igv hide variants that are missing IGV snapshot (default: FALSE)
#' @param layout slide layout choice: "individual" or "multiple" designed for a single or multiple individuals (default: "individual")
create_candidate_slides_ppt <- function(candidates, output_dir, output_cols,
                                        genemap2=NULL,
                                        GTEx_median_rpkm=NULL,
                                        GTEx_tissues=NULL,
                                        title_col = NULL,
                                        add_data_col = NULL,
                                        slide_template = NULL){
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
    panel_width <- 7.5
    omim_width <- 5
    
    candidates %>% 
        arrange(`inheritance model`, gene) %>% 
        split.data.frame(seq_len(nrow(.))) %>% 
        walk(function(cand) {
            gene <- cand$gene
            # Information table
            info_ft <-
                select(cand, all_of(output_cols)) %>% 
                t() %>% 
                as.data.frame() %>% 
                rownames_to_column('name') %>% 
                flextable() %>% 
                delete_part(part = "header") %>% 
                theme_zebra() %>% 
                italic(j = 1) %>% 
                bold(j = 1) %>% 
                colformat_char(j = 1, suffix = ':') %>% 
                align(j = 1, align = 'right', part = 'all') %>% 
                autofit() %>% 
                fit_to_width(seg_width)
            
            omim_ft <- NULL
            omim_table <- omim_table(gene, genemap2=genemap2, wrap = FALSE)
            if (!is.null(omim_table)) {
                omim_ft <-
                    omim_table %>% 
                    flextable() %>% 
                    theme_zebra() %>% 
                    autofit() %>% 
                    fit_to_width(omim_width)
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
                        fit_to_width(panel_width)
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
                ph_with(value = gtex_plot,
                        location = ph_location_template(left = seg_left[3],
                                                        top = top,
                                                        width = seg_width,
                                                        height = h1)) %>% 
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
                               location = ph_location_template(left = seg_left[1] + omim_width + 0.1,
                                                               top = top + h1 + pad_h,
                                                               width = panel_width))
                )} 
        })
    
    slides %>% print(target = file.path(output_dir, 'candidate_variants_report.pptx'))
}

