#!/bin/bash
# --------------------------------------------------------------------------------
# Example script to annotate VCF file using vcfanno and ANNOVAR
# --------------------------------------------------------------------------------

# Specify filenames
VCF_INPUT_FILE=/path/to/project/.../variants/variants.vcf.gz
OUTPUT_DIR=/path/to/project/.../annotated_variants/
OUTPUT_FILENAME_BASE=annotated_variants

# --------------------------------------------------------------------------------

# Specify location of vcfanno executable and configuration file
# (https://github.com/brentp/vcfanno)
VCFANNO_EXEC=/home/.../vcfanno/vcfanno_linux64
VCFANNO_TOML=vcfanno_conf_exome.toml

# Specify location of ANNOVAR script and humandb folder
# (http://annovar.openbioinformatics.org/en/latest/)
ANNOVAR_SCRIPT=/home/.../annovar/table_annovar.pl
ANNOVAR_HUMANDB_DIR=/home/.../annovar/humandb/

# --------------------------------------------------------------------------------

# Create directory and annotate VCF using vcfanno and ANNOVAR

# OUTPUT_DIR=$(dirname "${OUTPUT_FILENAME_BASE}")
mkdir $OUTPUT_DIR
cd $OUTPUT_DIR

$VCFANNO_EXEC -p 8 $VCFANNO_TOML $VCF_INPUT_FILE > $OUTPUT_FILENAME_BASE.vcfanno.vcf

bgzip $OUTPUT_FILENAME_BASE.vcfanno.vcf
tabix -p vcf $OUTPUT_FILENAME_BASE.vcfanno.vcf.gz

perl $ANNOVAR_SCRIPT $OUTPUT_FILENAME_BASE.vcfanno.vcf.gz \
     $ANNOVAR_HUMANDB_DIR -buildver hg19 \
     -vcfinput -out $OUTPUT_FILENAME_BASE -remove \
     -protocol refGene,exac03,gnomad_exome,gnomad_genome,avsnp147,dbnsfp33a -operation g,f,f,f,f,f -nastring .

bgzip $OUTPUT_FILENAME_BASE.hg19_multianno.vcf
tabix -p vcf $OUTPUT_FILENAME_BASE.hg19_multianno.vcf.gz
gzip $OUTPUT_FILENAME_BASE.avinput
gzip $OUTPUT_FILENAME_BASE.hg19_multianno.txt

