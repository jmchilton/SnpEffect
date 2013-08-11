#!/bin/sh

#-------------------------------------------------------------------------------
# Requirements:
#	- A computer (laptop, desktop). There is no need for big computers or clusters
#	- Operating system: Unix, Linux, Mac OS (or Windows + CygWin)
#	- Memory requirements: 4GB
#	- Java version 1.6 or higher (most modern computers have it)
#	- Hard drive space: 1GB
#
# Note:
#
# The sample dataaset consists of finding the generic cause of Cystic fibrosis.
# Given that this is an excercise and due to restricions in using 
# human sequencing data, we constructed the dataset as follows:
#	i) We used real sequencing data from "CEPH_1463" dataset. The dataset
#	   is provided by Complete Genomics diversity panel. It consists of
#	   sequencing of a family: 4 grandparents, 2 parents and 11 siblings
#
#	ii) In that dataset, we added a know Mendelian dicease mutation on three siblings (cases)
#
#	iii) We made this change while being consistent with the underlying heplotype structure.
#
# Goal:
#	The goal is to find the mutation causeing the mendelian rececive 
#   trait (Cystic Fibrosis).
#
# The assumptions are:
#	i) There is a mendelian recesive disorder: Three siblings are affected (cases), but parents and grand-parents are not.
#
#	ii) We assume this is caused by a high impact mutation in gene.
#
#-------------------------------------------------------------------------------

#---
# Download and install SnpEff and sample data
# Approximate time: 5 minutes
#---

# # Move to home directory
# cd

# # Download and install SnpEff
# curl -v -L http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip > snpEff_latest_core.zip
# unzip snpEff_latest_core.zip

# # Download and install SnpEff's Human Genome database
# java -jar snpEff.jar download -v GRCh37.71

# # Download sample data
# curl -v -L http://sourceforge.net/projects/snpeff/files/protocols_sample.zip > protocols_sample.zip
# unzip protocols_sample.zip

#---
# Step 1 : Annotate sample data
# Approximate time: 3 minutes
#---

# Annotate using human genome
# Notice that we use the '-lof' command line option to add loss of function and nonsense mediated decay tags

echo "Annotating VCF file"
java -Xmx4g -jar snpEff.jar \
  -v \
  -lof \
  -motif \
  -nextProt \
  GRCh37.71 \
  protocols_sample/chr7.vcf \
  > protocols_sample/chr7.eff.vcf

# Open HTML summary in browser, look for QC indicators

#---
# Step 2 : Count variants in case and control groups
# Approximate time: 1 hour
#---

# We know that there are three cases (NA12879, NA12885, NA12886). We created 
# a TFAM file accordingly. 
# Using "SnpSift caseControl" we count the number homozygous, heterozygous 
# and allele count for both cases and controls.

echo "Counting cases and controls"
java -Xmx1g -jar SnpSift.jar \
  caseControl \
  -v \
  -tfam protocols_sample/pedigree.tfam \
  protocols_sample/chr7.eff.vcf \
  > protocols_sample/chr7.eff.cc.vcf


# Note:  The program also calculates 
# p-values using different models. In this case, we know that the numner 
# of samples is not enough to reach genome wide significance. We will 
# not use those p-values in this excercise

#---
# Step 3 : Filter variants
# Approximate time: 1 minute
#---

# We will use "SnpSift filter" to reduce the number of candidate loci.
# The expression we use to filter is:
# 	i) We expect three cases to have a homozygous mutation. The INFO tag added 
#	in the previous step looks like this "Cases=1,1,3" where the first numbers 
#	represent the number of homozygous, heterozygous and allele count
#	So we want the first number "Cases[0]" (i.e. number of homozygous samples) 
#	to be 3:
#
#		(Cases[0] = 3)
#
#	ii) We expect none of the controls to have a homozygous mutation. The numbers
#	are represented as before (e.g. "Controls=8,6,22"). Since we want the first 
#	number to be zero, we write:
#
#		(Controls[0] = 0)
#
#	iii) Finally, we expect this to ba a high impact mutation. Since we know that
#	there are several effect per variant, we ask for ANY effect (EFF[*]) to have an
#	impact 'HIGH'. So the filter expression is
#
#		(EFF[*].IMPACT = 'HIGH')
#
# Putting all this toghether, we get:

echo "Filtering variants"
cat protocols_sample/chr7.eff.cc.vcf \
  | java -jar SnpSift.jar filter \
    "(Cases[0] = 3) & (Controls[0] = 0) & (EFF[*].IMPACT = 'HIGH')" \
  > protocols_sample/chr7.filtered.hom.high.vcf

echo "Show variants in a simple format"
cat protocols_sample/chr7.filtered.hom.high.vcf \
	| ./scripts/vcfInfoOnePerLine.pl \
	| less -S

# Only one result, show effects

# Plot pedigree
echo "Plot variant in pedigree"
java -jar SnpSift.jar pedShow \
	protocols_sample/pedigree.tfam \
	protocols_sample/chr7.filtered.hom.high.vcf \
	protocols_sample/chart

# Open protocols_sample/chart/7_117227832/index.html in your browser

#---
# Optional step: Annotating for GATK pipelines
#
# This is an optional step showing how to include 
# annotations into GATK's pipeline
#
# References:
#	GATK:							http://www.broadinstitute.org/gatk/
#   SnpEff annotations in GATK:		http://snpeff.sourceforge.net/SnpEff_manual.html#gatk
#---

# WARNING: These variables have to be changed according
#          where programs and files are installed in your 
#          system
ref=$HOME/snpEff/data/genomes/GRCh37.71.fa
dict=$HOME/snpEff/data/genomes/GRCh37.71.dict

gatk=$HOME/tools/gatk/GenomeAnalysisTK.jar 
picard=$HOME/tools/picard/

# Annotate using SnpEff for GATK
echo "Annotating VCF file for GATK"
java -Xmx4g -jar snpEff.jar \
  -v \
  -o gatk \
  GRCh37.71 \
  protocols_sample/chr7.vcf \
  > protocols_sample/chr7.eff.gatk.vcf

# # Create genome index file (if needed)
# echo "Indexing Genome reference FASTA file: $ref"
# samtools faidx $ref
# 
# # Create dictionary file (if needed)
# echo "Creating Genome reference dictionary file: $dict"
# java -jar $picard/CreateSequenceDictionary.jar R= $ref O= $dict

# # Use GATK's VariantAnnotator to add SnpEff's annotations
echo "Using GATK's VariantAnnotator. Path to GATK: $gatk"
java -Xmx4g -jar $gatk \
    -T VariantAnnotator \
    -R $ref \
    -A SnpEff \
    --variant protocols_sample/chr7.vcf \
    --snpEffFile protocols_sample/chr7.eff.gatk.vcf \
    -L protocols_sample/chr7.vcf \
    -o protocols_sample/chr7.gatk.vcf

