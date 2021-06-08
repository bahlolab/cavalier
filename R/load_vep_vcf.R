#' Load variants from VCF file (ANN field from VEP, ExAC counts by vcfanno)
#' 
#' @param vcf_filename VEP annotated VCF filename
#' @param sampleID list of IDs and names for samples of interest
# #' @examples
# #' ***TODO***

# assume vcf does not contain multiallelic variants
# require PICK column annotation by VEP
load_vep_vcf <- function(vcf_filename, sampleID, field = 'CSQ') {
    vep_col_spec <-
        readr::cols(
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
            HGNC_ID = col_double(),
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
    
    vcf <- vcfR::read.vcfR(vcf_filename)
    vcf_tidy <- vcfR::vcfR2tidy(vcf)
    
    # Get ANN format column names
    ANN_columns <- strsplit(strsplit(vcf_tidy$meta$Description[vcf_tidy$meta$ID == field], "Format: ")[[1]][2], "|", fixed=TRUE)[[1]]
    
    # get vep annotation, index corresponds variant record index
    # filter for VEP picked annotation (i.e. CANONICAL)
    ANN <-
        tibble(index = seq_along(vcf_tidy$fix[[field]]),
               data = str_split(vcf_tidy$fix[[field]], ',')) %>% 
        unnest(data) %>% 
        (function(x) {
            x <- str_split(x$data, '\\|', simplify = TRUE) %>% 
                magrittr::set_colnames(ANN_columns) %>% 
                as_tibble() %>% 
                readr::type_convert(col_types = vep_col_spec) %>% 
                { bind_cols(select(x, index), .)}
        }) %>% 
        filter(PICK == 1)

    # Combine and clean VCF data
    vars <-
        vcf_tidy$fix %>% 
        select(c("CHROM", "POS", "REF", "ALT", "QUAL", "FILTER", "DP", "QD", "MQ", "FS", "SOR", "MQRankSum", "ReadPosRankSum", "InbreedingCoeff", "AC")) %>% 
        (function(x) { bind_cols(x[ANN$index, ], ANN) }) %>% 
        rename(chromosome = CHROM,
               position = POS,
               reference = REF,
               alternate = ALT,
               change = Consequence,
               gene = SYMBOL,
               dbSNP = Existing_variation,
               Polyphen2 = PolyPhen,
               MAF_1000G = AF,
               MAF_gnomAD = gnomAD_AF) %>% 
        mutate(MAF_1000G = replace_na(MAF_1000G, 0),
               MAF_gnomAD = replace_na(MAF_gnomAD, 0),
               SIFT_score = str_extract(SIFT, '(?<=\\()[0-9\\.]+(?=\\)$)') %>% as.numeric(),
               SIFT = str_extract(SIFT, '^.+(?=\\([0-9\\.]+\\)$)') %>% str_remove('_low_confidence'),
               Polyphen2_score = str_extract(Polyphen2, '(?<=\\()[0-9\\.]+(?=\\)$)') %>% as.numeric(),
               Polyphen2 = str_extract(Polyphen2, '^.+(?=\\([0-9\\.]+\\)$)')) %>% 
        select(., chromosome, position, reference, alternate, gene,
               starts_with('MAF'), starts_with('SIFT'), starts_with('Polyphen2'), everything())
    
    for (sID in names(sampleID)) {
        vcf_tidy_gt_sID <- vcf_tidy$gt[vcf_tidy$gt$Indiv == sID, ]
        sID_name <- sampleID[[sID]]
        vars <-
            vars %>% 
            mutate(!!paste(sID_name, "genotype") := vcf_tidy_gt_sID$gt_GT[index],
                   !!paste(sID_name, "GT quality") := vcf_tidy_gt_sID$gt_GQ[index],
                   !!paste(sID_name, "depth (R,A)") := vcf_tidy_gt_sID$gt_AD[index])
    }
    
    vars$chromosome <- ifelse(startsWith(vars$chromosome, "chr"), vars$chromosome, paste0("chr", vars$chromosome))
    
    # Construct end position from start position and difference between length of ref and alt alleles
    length_ref <- sapply(vars$reference, nchar)
    length_alt <- sapply(vars$alternate, nchar)
    # Start and end positions are the same for insertions or SNV but different for deletions
    vars$end <- ifelse(length_ref <= length_alt, vars$position, vars$position + length_ref - length_alt - 1)
    
    vars$change <- gsub("&", ";", vars$change, fixed=TRUE)
    
    # Annotation  *** TODO: add
    AA_split <- strsplit(ANN$Amino_acids, "/", fixed=TRUE)
    AA1 <- sapply(AA_split, function(x){x[1]})
    AA2 <- sapply(AA_split, function(x){x[2]})
    AA1[is.na(AA1)] <- ""
    AA2[is.na(AA2)] <- ""
    vars$annotation <- paste0("p.", AA1, vcf_tidy$ANN$Protein_position, AA2)
    vars$annotation <- gsub(",", "", vars$annotation, fixed=TRUE)
    vars$annotation[vars$annotation == "p."] <- ""
    
    vars$dbSNP <- sapply(strsplit(vars$dbSNP, "&"), function(x){paste(x[startsWith(x, "rs")], collapse=";")})
    # fix gene names
    genes_wrong <- which(vars$gene %in% HGNC_alias$alias)
    vars$gene[genes_wrong] <- HGNC_alias$symbol[match(vars$gene[genes_wrong], HGNC_alias$alias)]
    # add Intolerance scores
    vars$RVIS <- rvis_exac_percentile(vars$gene)
    vars$GeVIR <- gevir_percentile(vars$gene)
    vars$LOEUF <- loeuf_percentile(vars$gene)
    vars$Grantham <- grantham_score(vars$annotation)
    
    return(vars)
}


