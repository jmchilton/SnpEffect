#!/bin/sh

VERSION_SNPEFF=2.2
VERSION_SNPSIFT=1.7

#---
# Build SnpEff
#---

cd $HOME/workspace/SnpEff/

mvn assembly:assembly
cp target/snpEff-$VERSION_SNPEFF-jar-with-dependencies.jar $HOME/snpEff/snpEff.jar


# Install JAR file in local Maven repo
mvn install:install-file \
	-Dfile=target/snpEff-$VERSION_SNPEFF-jar-with-dependencies.jar \
	-DgroupId=ca.mcgill.mcb.pcingola \
	-DartifactId=snpEff \
	-Dversion=2.2 \
	-Dpackaging=jar

cd - 

#---
# Build SnpSift
#---
cd $HOME/workspace/SnpSift/

mvn assembly:assembly
cp target/snpSift-$VERSION_SNPSIFT-jar-with-dependencies.jar $HOME/snpEff/SnpSift.jar

cd - 

