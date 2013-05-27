#!/bin/sh -e

source `dirname $0`/config.sh

#mkdir download
cd download

# #---
# # Download from ENSEMBL
# #---
# 
# # Download GTF files (annotations)
# wget -r -A "*gtf.gz" "ftp://ftp.ensembl.org/pub/release-$ENSEMBL_RELEASE/gtf/"
# 
# # Download FASTA files (reference genomes)
# wget -r -A "*dna.toplevel.fa.gz" "ftp://ftp.ensembl.org/pub/release-$ENSEMBL_RELEASE/fasta/"
# 
# # Download CDS sequences
# wget -r -A "*cdna.all.fa.gz" "ftp://ftp.ensembl.org/pub/release-$ENSEMBL_RELEASE/fasta/"
# 
# # Download PROTEIN sequences
# wget -r -A "*.pep.all.fa.gz" "ftp://ftp.ensembl.org/pub/release-$ENSEMBL_RELEASE/fasta/"
# 
# # Download regulation tracks
# wget -r -A "*AnnotatedFeatures.gff.gz" "ftp://ftp.ensembl.org/pub/release-$ENSEMBL_RELEASE/regulation/"
# wget -r -A "*MotifFeatures.gff.gz" "ftp://ftp.ensembl.org/pub/release-$ENSEMBL_RELEASE/regulation/"
# 
# #---
# # Create directory structure
# #---
# 
# # Move all GTF and FASTA downloaded files to this directory
# mv `find ftp.ensembl.org -type f -iname "*gtf.gz" -or -iname "*.fa.gz"` .
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
# mkdir -p data/GRCh37.$ENSEMBL_RELEASE/
# cp ftp.ensembl.org/pub/release-$ENSEMBL_RELEASE/regulation/homo_sapiens/AnnotatedFeatures.gff.gz data/GRCh37.$ENSEMBL_RELEASE/regulation.gff.gz
# cp ftp.ensembl.org/pub/release-$ENSEMBL_RELEASE/regulation/homo_sapiens/MotifFeatures.gff.gz data/GRCh37.$ENSEMBL_RELEASE/motif.gff.gz
# 
# mkdir -p data/GRCm38.$ENSEMBL_RELEASE/
# cp ftp.ensembl.org/pub/release-$ENSEMBL_RELEASE/regulation/mus_musculus/AnnotatedFeatures.gff.gz data/GRCm38.$ENSEMBL_RELEASE/regulation.gff.gz
# cp ftp.ensembl.org/pub/release-$ENSEMBL_RELEASE/regulation/mus_musculus/MotifFeatures.gff.gz data/GRCm38.$ENSEMBL_RELEASE/motif.gff.gz
# 
# #---
# # Config file entries
# #---
# 
# (
# for fasta in *.cdna.all.fa.gz
# do
# 	genome=`../scripts/file2GenomeName.pl $fasta | cut -f 4`
# 	short=`../scripts/file2GenomeName.pl $fasta | cut -f 5`
# 
# 	# Individual genome entry
# 	echo -e "$short.genome : $genome"
# 	echo -e "$short.reference : ftp://ftp.ensembl.org/pub/release-$ENSEMBL_RELEASE/gtf/"
# 	echo
# done
# ) | tee ../snpEff.ensembl.$ENSEMBL_RELEASE.config
# 
# #---
# # ENSEMBL is compressing files in a way that is not fully compatible with Java's gzip library
# #---
# 
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
# 
# #---
# # Back to parent dir
# #---
# cd - > /dev/null
# 
# #---
# # Move data dir to 'real' data dir
# #---
# 
# mv download/data/genomes/* data/genomes/
# rmdir download/data/genomes
# 
# mv download/data/* data/
# 
# echo "Done!"
# 
