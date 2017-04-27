### Example cavalier script

## Specify Project Information and Filter Settings
sampleID <- list("XXXXX"="proband")
pedigree <- "/path/to/family.ped"
ANNOVAR_vcf <- "/path/to/annovar.vcf.gz"
output_dir <- "/path/to/cavalier/output/"

# List of inheritance models and corresponding MAF thresholds for each
inheritance_MAF <- list("individual dominant"  = 0.0001,
                        "individual recessive" = 0.01,
                        "individual comp het"  = 0.01)

# Need to register and download GTEx tissue median RPKM table
# GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct.gz
# from https://gtexportal.org/home/datasets and specify filename for GTEx function
GTEx_median_rpkm_file <- "/path/to/GTEx_data/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct.gz"

# Need to register and download OMIM genemap2.txt table from 
# https://omim.org/downloads/ and specify filename here for OMIM function
OMIM_genemap2_file <- "/path/to/OMIM_data/genemap2.txt"


# Use devtools to install cavalier R package
library(devtools)
install_github("bahlolab/cavalier")

library(cavalier)

## Check pedigree relatedness, sex, ancestry with peddy (optional step; requires peddy package)
## [https://github.com/brentp/peddy]
run_peddy(output_dir, pedigree, ANNOVAR_vcf)

## Load variants (currently different functions for VCF files produced by VEP or ANNOVAR)
#  [Needs to be tested more broadly to avoid errors with different annotation settings]
vars <- load_annovar_vcf(ANNOVAR_vcf, sampleID)
# vars <- load_vep_vcf(vcf_filename, sampleID)

## Filter variants for quality
qualvars <- quality_filter_variants(vars)

# Filter variants based on inheritance
# (see function for additional options: include/exclude particular genes, regions, or types of change, etc)
candvars <- filter_variants(qualvars, inheritance_MAF, sampleID, 
                            region_include=c("exonic", "splicing", "exonic;splicing"),
                            change_exclude=c("synonymous SNV"))

## Create IGV batch script [output_dir]/data/igv_batch.txt
candvars <- create_igv_batch_script(candvars, output_dir)

# *** NEED TO MANUALLY LOAD IGV AND RUN IGV BATCH SCRIPT TO CREATE IGV SNAPSHOTS ***

## Create cavalier output
create_cavalier_output(candvars, output_dir, sampleID, hide_missing_igv=TRUE, layout="individual", 
                       genemap2=OMIM_genemap2_file, GTEx_median_rpkm=GTEx_median_rpkm_file)

