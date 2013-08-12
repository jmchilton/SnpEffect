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
	-interval protocols/ex2_regulatory.bed \
	GRCh37.71 \
	protocols/ex2.vcf \
	> protocols/ex2.eff.vcf

#---
# Step 2: Add conservation scores
#---

java -Xmx1g -jar SnpSift.jar \
	phastCons \
	-v \
	protocols/phastcons \
	protocols/ex2.eff.vcf \
	> protocols/ex2.eff.cons.vcf

#---
# Step 3: Filter variants
#---

cat protocols/ex2.eff.cons.closest.vcf \
  | java -jar SnpSift.jar filter \
    "(EFF[*].EFFECT = 'CUSTOM[ex2_regulatory]') & (exists PhastCons) & (PhastCons > 0.9)" \
  > protocols/ex2.filtered.vcf

