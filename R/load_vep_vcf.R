#' @importFrom stringr str_c str_extract str_split_fixed str_detect
#' @importFrom magrittr '%>%' set_colnames
#' @importFrom dplyr transmute coalesce
get_vep_ann <- function(gds, vep_field,
                        add_annot = character()) {
    
    vep_ann_names <- 
        SeqArray::header(gds)$INFO %>%
        as.data.frame() %>% 
        tibble::rownames_to_column() %>% 
        as_tibble() %>%
        filter(str_detect(rowname, vep_field),
               str_detect(Description, 'Ensembl VEP'))%>%
        pull(Description) %>%
        str_remove('.+Format: ') %>% 
        str_split('\\|', simplify = T) %>% 
        c()
    
    vep_ann <- seqGetData(gds, str_c('annotation/info/', vep_field))
    vid <- tibble(variant_id = rep(seqGetData(gds, 'variant.id'), times = vep_ann$length))
    
    vep_raw <-
        str_split_fixed(vep_ann$data, '\\|', length(vep_ann_names)) %>%
        set_colnames(vep_ann_names) %>% 
        as_tibble() %>% 
        readr::type_convert(col_types = vep_col_spec)
    
    vep_clean <-
        vep_raw %>% 
        transmute(gene = coalesce(hgnc_ensembl2sym(Gene),
                                  hgnc_id2sym(HGNC_ID),
                                  hgnc_sym2sym(SYMBOL)),
                  hgnc_id = HGNC_ID,
                  ensembl_gene = Gene,
                  ensembl_transcript = str_extract(Feature, '^ENST.+'),
                  ensembl_protein = str_extract(ENSP, '^ENSP.+'),
                  consequence = Consequence,
                  impact = ordered(IMPACT, c('MODIFIER', 'LOW', 'MODERATE', 'HIGH')),
                  hgvs_genomic = HGVSg,
                  hgvs_coding = str_extract(HGVSc, '(?<=:).+$'),
                  hgvs_protein = str_extract(HGVSp, '(?<=:).+$') %>% str_replace('%3D', '='),
                  db_snp = str_extract(Existing_variation, 'rs[0-9]+'),
                  af_gnomad = gnomAD_AF,
                  af_1000G = AF,
                  af_popmax = MAX_AF,
                  sift = str_extract(SIFT, '^.+(?=\\([0-9\\.]+\\)$)') %>% str_remove('_low_confidence'),
                  sift_score = str_extract(SIFT, '(?<=\\()[0-9\\.]+(?=\\)$)') %>% as.numeric(),
                  polyphen = str_extract(PolyPhen, '^.+(?=\\([0-9\\.]+\\)$)'),
                  polyphen_score = str_extract(PolyPhen, '(?<=\\()[0-9\\.]+(?=\\)$)') %>% as.numeric()) %T>% 
        with(assert_that(all(sift %in% c(NA, 'tolerated', 'deleterious'))),
             assert_that(all(polyphen %in% c(NA, 'benign', 'possibly_damaging', 'probably_damaging', 'unknown')))) %>% 
        mutate(sift = ordered(sift,  c('tolerated', 'deleterious')),
               polyphen = ordered(polyphen, c('benign', 'possibly_damaging', 'probably_damaging'))) %>% 
        # add rvis_percentile
        (function(data) `if`('rvis_percentile' %in% add_annot,
                             mutate(data, rvis_percentile = sym2rvis(gene)),
                             data)) %>% 
        # add gevir_percentile
        (function(data) `if`('gevir_percentile' %in% add_annot,
                             mutate(data, gevir_percentile = coalesce(sym2gevir(gene), ensembl2gevir(ensembl_gene))),
                             data)) %>% 
        # add loeuf_percentile
        (function(data) `if`('loeuf_percentile' %in% add_annot,
                             mutate(data, loeuf_percentile = coalesce(sym2loeuf(gene), ensembl2loeuf(ensembl_gene))),
                             data)) %>% 
        # add grantham_score
        (function(data) `if`('grantham_score' %in% add_annot,
                             mutate(data, grantham_score = grantham_score(hgvs_protein)),
                             data))
    
    return( bind_cols(vid, vep_clean) ) 
}

#' @importFrom readr cols col_character col_double col_integer
vep_col_spec <- cols(
    Allele = col_character(),
    Consequence = col_character(),
    IMPACT = col_character(),
    SYMBOL = col_character(),
    Gene = col_character(),
    Feature_type = col_character(),
    Feature = col_character(),
    BIOTYPE = col_character(),
    EXON = col_character(),
    INTRON = col_character(),
    HGVSc = col_character(),
    HGVSp = col_character(),
    cDNA_position = col_character(),
    CDS_position = col_character(),
    Protein_position = col_character(),
    Amino_acids = col_character(),
    Codons = col_character(),
    Existing_variation = col_character(),
    ALLELE_NUM = col_double(),
    DISTANCE = col_character(),
    STRAND = col_double(),
    FLAGS = col_character(),
    VARIANT_CLASS = col_character(),
    SYMBOL_SOURCE = col_character(),
    HGNC_ID = col_character(),
    CANONICAL = col_character(),
    MANE_SELECT = col_character(),
    MANE_PLUS_CLINICAL = col_character(),
    TSL = col_character(),
    APPRIS = col_character(),
    CCDS = col_character(),
    ENSP = col_character(),
    SWISSPROT = col_character(),
    TREMBL = col_character(),
    UNIPARC = col_character(),
    UNIPROT_ISOFORM = col_character(),
    GENE_PHENO = col_double(),
    SIFT = col_character(),
    PolyPhen = col_character(),
    DOMAINS = col_character(),
    miRNA = col_character(),
    AF = col_double(),
    AFR_AF = col_double(),
    AMR_AF = col_double(),
    EAS_AF = col_double(),
    EUR_AF = col_double(),
    SAS_AF = col_double(),
    AA_AF = col_double(),
    EA_AF = col_double(),
    gnomAD_AF = col_double(),
    gnomAD_AFR_AF = col_double(),
    gnomAD_AMR_AF = col_double(),
    gnomAD_ASJ_AF = col_double(),
    gnomAD_EAS_AF = col_double(),
    gnomAD_FIN_AF = col_double(),
    gnomAD_NFE_AF = col_double(),
    gnomAD_OTH_AF = col_double(),
    gnomAD_SAS_AF = col_double(),
    MAX_AF = col_double(),
    MAX_AF_POPS = col_character(),
    CLIN_SIG = col_character(),
    SOMATIC = col_character(),
    PHENO = col_character(),
    PUBMED = col_character(),
    CHECK_REF = col_character(),
    MOTIF_NAME = col_character(),
    MOTIF_POS = col_double(),
    HIGH_INF_POS = col_character(),
    MOTIF_SCORE_CHANGE = col_character(),
    TRANSCRIPTION_FACTORS = col_character(),
    PICK = col_integer()
)