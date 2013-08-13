#!/bin/sh -e

# # Update ClinVar
# cd db/clinVar
# wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf/clinvar_00-latest.vcf.gz
# gunzip clinvar_00-latest.vcf.gz 
# cd -

# # Update GwasCatalog
# cd db/gwasCatalog
# wget http://www.genome.gov/admin/gwascatalog.txt
# cd -

# # Update dbSnp
# cd db/dbSnp/
# wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz
# gunzip 00-All.vcf.gz
# mv 00-All.vcf dbSnp.vcf
# cd -

#---
# Create a big ZIP file
#---
ZIP=db_GRCh.`date +"%Y%m%d"`.zip
zip -r $ZIP \
    db/dbSnp/ \
    db/encode \
    db/epigenome/Histone \
    db/gwasCatalog \
    db/jaspar \
    db/miRNA \
    db/msigDb/*.gmt \
    db/phastCons/hg19/ \

# Upload to sourceforge
scp $ZIP pcingola,snpeff@frs.sourceforge.net:/home/frs/project/s/sn/snpeff/databases/
