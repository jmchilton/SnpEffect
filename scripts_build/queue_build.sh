#!/bin/sh

#---
# Create build queue entries
#---
echo Creating build queue 'queue_build.txt' file
rm -vf queue_build.txt

# Build from refseq files
echo "java -Xmx10G -jar snpEff.jar build -v -refseq hg19" >> queue_build.txt

# Build from TXT files
for genes in `ls data/*/genes.txt* | grep -v hg19`
do
	dir=`dirname $genes`
	genomeName=`basename $dir`
	echo "java -Xmx10G -jar snpEff.jar build -v -txt $genomeName"
done | sort >> queue_build.txt

# Build from GFF2 files
echo "java -Xmx10G -jar snpEff.jar build -v -gff2 amel2" >> queue_build.txt

# Build from GFF3 files
for genes in `ls data/*/genes.gff* | grep -v amel2`
do
	dir=`dirname $genes`
	genomeName=`basename $dir`
	echo "java -Xmx10G -jar snpEff.jar build -v -gff3 $genomeName"
done | sort >> queue_build.txt

# Build from GTF22 files
for genes in data/*/genes.gtf*
do
	dir=`dirname $genes`
	genomeName=`basename $dir`
	echo "java -Xmx10G -jar snpEff.jar build -v -gtf22 $genomeName"
done | sort >> queue_build.txt

# Build from GenBank files
for genes in data/*/genes.gb*
do
	dir=`dirname $genes`
	genomeName=`basename $dir`
	echo "java -Xmx10G -jar snpEff.jar build -v -genbank $genomeName"
done | sort >> queue_build.txt

#---
# Build ALL genomes
#---
./scripts/queue.pl 24 24 1 queue_build.txt

#---
# Create build reports
#---
echo "Creating build report : cds_check.txt"
grep "CDS check" *.stdout | cut -f 3- -d : | grep -v "^$" | sort | tee cds_check.txt

echo "Creating build report : protein_check.txt"
grep "Protein check" *.stdout | cut -f 3- -d : | grep -v "^$" | sort | tee protein_check.txt

