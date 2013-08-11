#!/bin/sh

#-------------------------------------------------------------------------------
# Requirements:
#	- A computer (laptop, desktop). There is no need for big computers or clusters
#	- Operating system: Unix, Linux, Mac OS (or Windows + CygWin)
#	- Memory requirements: 4GB
#	- Java version 1.6 or higher (most modern computers have it)
#	- Hard drive space: 1GB
#-------------------------------------------------------------------------------

#---
# Step 1 : Annotate sample data
# Approximate time: 3 minutes
#---

# Annotate using human genome
# Notice that we use the '-motif' command line option to add TBFS information

java -Xmx4g -jar snpEff.jar \
	-motif \
	-interval protocols_sample/tbx5_regulatory.bed \
	GRCh37.71 \
	protocols_sample/tbx5.vcf \
	> protocols_sample/tbx5.eff.vcf

#---
# Step 2: Add conservation scores
#---

java -Xmx1g -jar SnpSift.jar \
	phastCons \
	-v \
	protocols_sample/phastcons \
	protocols_sample/tbx5.eff.vcf \
	> protocols_sample/tbx5.eff.cons.vcf

#---
# Step 3: Filter variants
#---

cat protocols_sample/tbx5.eff.cons.closest.vcf \
  | java -jar SnpSift.jar filter \
    "(EFF[*].EFFECT = 'CUSTOM[tbx5_regulatory]') & (exists PhastCons) & (PhastCons > 0.9)" \
  > protocols_sample/tbx5.filtered.vcf

