#!/bin/sh

#-------------------------------------------------------------------------------
#
# Upload files to SourceForge web 
#
#																Pablo Cingolani
#-------------------------------------------------------------------------------

# Include variables
source `dirname $0`/config.sh

#---
# Upload to ZIP files and databases
#---
# Core program
scp snpEff_v${SNPEFF_VERSION}_core.zip pcingola,snpeff@frs.sourceforge.net:/home/frs/project/s/sn/snpeff/snpEff_latest_core.zip
scp snpEff_v${SNPEFF_VERSION}_core.zip pcingola,snpeff@frs.sourceforge.net:/home/frs/project/s/sn/snpeff/
		
# Individual databases
scp snpEff_v${SNPEFF_VERSION}_*.zip pcingola,snpeff@frs.sourceforge.net:/home/frs/project/s/sn/snpeff/databases/v${SNPEFF_VERSION}/

# SnpSift
scp SnpSift.jar pcingola,snpeff@frs.sourceforge.net:/home/frs/project/s/sn/snpeff/SnpSift_v${SNPSIFT_VERSION}.jar
scp SnpSift.jar pcingola,snpeff@frs.sourceforge.net:/home/frs/project/s/sn/snpeff/SnpSift_latest.jar

#---
# Update SnpEff web pages
#---

# Create versions file (html/versions.txt)
./scripts_build/versions.sh

# Upload HTML pages
cd $HOME/workspace/SnpEff/html/

# Copy html and txt files 
scp style.css *.html *.txt pcingola,snpeff@frs.sourceforge.net:htdocs/
		
# Copy images
scp -r  images/ pcingola,snpeff@frs.sourceforge.net:htdocs/images/

