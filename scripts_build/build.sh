#!/bin/sh

#-------------------------------------------------------------------------------
#
# SnpEff release process
#
#																Pablo Cingolani
#-------------------------------------------------------------------------------

# Include variables
source `dirname $0`/config.sh

# Create JAR files
./scripts_build/make.sh 

#---
# Download databases
#---

# ENSEMBL databases: 
./scripts_build/download_ensembl.sh
		
# hg19 database
./scripts_build/download_hg19.sh

# Download NCBI genomes
./scripts_build/download_ncbi.sh

#---
# Build databases
#---

# Delete old databases
rm -vf data/*/*.bin

# Create queue entries
./scripts_build/queue_build.sh

# Build 
./scripts/queue.pl 22 22 1 queue_build.txt

#---
# Build distribution files
#---
./scripts_build/distro.sh

#---
# Upload to SourceForge
#---
./scripts_build/uploadSourceForge.sh

