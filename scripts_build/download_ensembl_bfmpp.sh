#!/bin/sh -e

#-------------------------------------------------------------------------------
# Download ENSEMBL's Bacteria, Fungi, Metazoa, Plants, Protist
#
# Note: This uses the alternative site ftp://ftp.ensemblgenomes.org
#       instead of the usual ftp://ftp.ensembl.org
#
#
#																Pablo Cingolani
#-------------------------------------------------------------------------------
source `dirname $0`/config.sh

#mkdir download
cd download

site="ftp://ftp.ensemblgenomes.org"

wget_wait=1
wget="wget --wait=$wget_wait -r -nc "

#---
# Download from ENSEMBL
#---

for org in bacteria fungi metazoa misc_data plants protists
do
	# Download GTF files (annotations)
	$wget -A "*gtf.gz" "$site/pub/$org/release-$ENSEMBL_BFMPP_RELEASE/gtf/"
	 
	# Download FASTA files (reference genomes)
	$wget -A "*dna.toplevel.fa.gz" "$site/pub/$org/release-$ENSEMBL_BFMPP_RELEASE/fasta/"

	# Download CDS sequences
	$wget -A "*cdna.all.fa.gz" "$site/pub/$org/release-$ENSEMBL_BFMPP_RELEASE/fasta/"

	# Download PROTEIN sequences
	$wget -A "*.pep.all.fa.gz" "$site/pub/$org/release-$ENSEMBL_BFMPP_RELEASE/fasta/"

done

# Download regulation tracks
#wget -r -A "*AnnotatedFeatures.gff.gz" "ftp://ftp.ensembl.org/pub/release-$ENSEMBL_BFMPP_RELEASE/regulation/"
#wget -r -A "*MotifFeatures.gff.gz" "ftp://ftp.ensembl.org/pub/release-$ENSEMBL_BFMPP_RELEASE/regulation/"

#---
# Create directory structure
#---

# # Move all downloaded file to this directory
# mv `find ftp.ensembl.org -type f` .
# 
# # Gene annotations files
# for gtf in *.gtf.gz
# do
# 	short=`../scripts/file2GenomeName.pl $gtf | cut -f 5`
# 	echo ANNOTATIONS: $short
# 
# 	mkdir -p data/$short
# 	cp $gtf data/$short/genes.gtf.gz
# done
#  
# # Reference genomes files
# mkdir -p data/genomes
# for fasta in *.dna.toplevel.fa.gz
# do
# 	genome=`../scripts/file2GenomeName.pl $fasta | cut -f 5`
# 	echo REFERENCE: $genome
# 
# 	cp $fasta data/genomes/$genome.fa.gz
# done
# 
# # CDS genomes files
# for fasta in *.cdna.all.fa.gz
# do
# 	genome=`../scripts/file2GenomeName.pl $fasta | cut -f 5`
# 	echo CDS: $genome
# 
# 	cp $fasta data/$genome/cds.fa.gz
# done
# 
# # Protein seuqence files
# for pep in *.pep.all.fa.gz
# do
# 	short=`../scripts/file2GenomeName.pl $pep | cut -f 5`
# 	echo PROTEIN: $short
# 
# 	mkdir -p data/$short
# 	cp $pep data/$short/protein.fa.gz
# done
#
# # Regunation tracks
# mkdir -p data/GRCh37.$ENSEMBL_BFMPP_RELEASE/
# mv ftp.ensembl.org/pub/release-$ENSEMBL_BFMPP_RELEASE/regulation/homo_sapiens/AnnotatedFeatures.gff.gz data/GRCh37.$ENSEMBL_BFMPP_RELEASE/regulation.gff.gz
# mv ftp.ensembl.org/pub/release-$ENSEMBL_BFMPP_RELEASE/regulation/homo_sapiens/MotifFeatures.gff.gz data/GRCh37.$ENSEMBL_BFMPP_RELEASE/motif.gff.gz
# 
# mkdir -p data/GRCm38.$ENSEMBL_BFMPP_RELEASE/
# mv ftp.ensembl.org/pub/release-$ENSEMBL_BFMPP_RELEASE/regulation/mus_musculus/AnnotatedFeatures.gff.gz data/GRCm38.$ENSEMBL_BFMPP_RELEASE/regulation.gff.gz
# mv ftp.ensembl.org/pub/release-$ENSEMBL_BFMPP_RELEASE/regulation/mus_musculus/MotifFeatures.gff.gz data/GRCm38.$ENSEMBL_BFMPP_RELEASE/motif.gff.gz

#---
# Config file entries
#---

# for fasta in *.cdna.all.fa.gz
# do
# 	genome=`../scripts/file2GenomeName.pl $fasta | cut -f 4`
# 	short=`../scripts/file2GenomeName.pl $fasta | cut -f 5`
# 
# 	# Individual genome entry
# 	echo -e "$short.genome : $genome"
# 	echo -e "$short.reference : ftp://ftp.ensembl.org/pub/release-$ENSEMBL_BFMPP_RELEASE/gtf/"
# 	echo
# done

#---
# ENSEMBL is files in a way that is not fully compatible with Java's gzip library
#---

# rm -vf queue_gunzip.txt queue_gunzip.txt
# for g in `find . -iname "*.gz"`
# do
# 	f=`dirname $g`/`basename $g .gz`
# 	echo "Un-compress / compress: $g"
# 	echo "gunzip -v $g" >> queue_gunzip.txt
# 	echo "gzip -v $f" >> queue_gzip.txt
# done
# 
# # Uncompress all files
# ../scripts/queue.pl 22 22 1 queue_gunzip.txt
# 
# # Compress all files
# ../scripts/queue.pl 22 22 1 queue_gzip.txt

#---
# Back to parent dir
#---
cd - > /dev/null

#---
# Move data dir to 'real' data dir
#---
# mv download/data/genomes/* data/genomes/
# mv download/data/* data/

