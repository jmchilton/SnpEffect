#!/bin/sh

#-------------------------------------------------------------------------------
#
# Create databases
#
#-------------------------------------------------------------------------------

#---
# Create build queue 
#---
snpeff="java -Xmx10G -jar snpEff.jar "
snpeffBuild="$snpeff build -v"

echo Creating build queue 'queue_build.txt' file
rm -vf queue_build.txt

# RefSeq files
echo "$snpeffBuild -refseq hg19" >> queue_build.txt

# TXT files
for genes in `ls data/*/genes.txt* | grep -v hg19`
do
	dir=`dirname $genes`
	genomeName=`basename $dir`
	echo "$snpeffBuild -txt $genomeName"
done | sort >> queue_build.txt

# GFF2 files
echo "$snpeffBuild -gff2 amel2" >> queue_build.txt

# GFF3 files
for genes in `ls data/*/genes.gff* | grep -v amel2`
do
	dir=`dirname $genes`
	genomeName=`basename $dir`
	echo "$snpeffBuild -gff3 $genomeName"
done | sort >> queue_build.txt

# GTF22 files
for genes in data/*/genes.gtf*
do
	dir=`dirname $genes`
	genomeName=`basename $dir`
	echo "$snpeffBuild -gtf22 $genomeName"
done | sort >> queue_build.txt

# GenBank files
for genes in data/*/genes.gb*
do
	dir=`dirname $genes`
	genomeName=`basename $dir`
	echo "$snpeffBuild -genbank $genomeName"
done | sort >> queue_build.txt

#---
# Build genomes
#---
./scripts/queue.pl 24 24 1 queue_build.txt

#---
# Special builds: Human
#---
$snpeff buildNextProt -v GRCh37.$ENSEMBL_RELEASE db/nextProt/

#---
# Create build reports
#---
echo "Creating build report : cds_check.txt"
grep "CDS check" *.stdout | cut -f 3- -d : | grep -v "^$" | sort | tee cds_check.txt

echo "Creating build report : protein_check.txt"
grep "Protein check" *.stdout | cut -f 3- -d : | grep -v "^$" | sort | tee protein_check.txt

