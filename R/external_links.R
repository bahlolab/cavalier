
gnomad_link <- function(variants,
                        ref_genome = get_cavalier_opt('ref_genome'))
{
  assert_that(is.data.frame(variants),
              is_scalar_character(ref_genome),
              ref_genome %in% c('hg38', 'hg19'))
  with(variants,
       assert_that(is.character(chrom),
                   is.integer(pos),
                   is.character(ref),
                   is.character(alt)))
  
  with(variants,
       str_c('https://gnomad.broadinstitute.org/variant/',
             str_c(str_remove(chrom, 'chr'), pos, ref, alt, sep = '-'),
             '?dataset=gnomad_',
             if_else(ref_genome == 'hg38', 'r3', 'r2_1')))
}

dbsnp_link <- function(rsid) 
{
  assert_that(is.character(rsid))
  
  str_c('https://www.ncbi.nlm.nih.gov/snp/',
        str_extract(rsid, '^rs[0-9]+$'))
}

genecards_link <- function(gene) 
{
  assert_that(is.character(gene))
  
  str_c('https://www.genecards.org/cgi-bin/carddisp.pl?gene=',
        gene)
}

ensembl_gene_link <- function(ensembl_gene) 
{
  assert_that(is.character(ensembl_gene))
  
  str_c('https://ensembl.org/Homo_sapiens/Gene/Summary?g=',
        str_extract(ensembl_gene, '^ENSG[0-9]+$'))
}

