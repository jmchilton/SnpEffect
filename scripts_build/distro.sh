#!/bin/sh

#------------------------------------------------------------------------------
# Create a zip file for distribution
# Note: Only binary data is included (no raw gene info / genomes)
#
#                                      Pablo Cingolani 2010
#------------------------------------------------------------------------------

source `dirname $0`/config.sh

DIR=snpEff_$SNPEFF_VERSION_REV
rm -rvf $DIR
mkdir $DIR

# Copy core files
cp snpEff.config snpEff.jar SnpSift.jar demo.1kg.vcf $DIR
# cp -rvfH galaxy scripts $DIR
cp -rvf galaxy scripts $DIR

cd $DIR
rm -rvf `find . -name "CVS" -type d`
cd -

# Change name to 'snpEff' (so that config file can be used out of the box)
mv $DIR snpEff

# Create 'core' zip file
ZIP="snpEff_v"$SNPEFF_VERSION_REV"_core.zip"
rm -f $ZIP 2> /dev/null
zip -r $ZIP snpEff

# Create ZIP file for each database
for d in `ls data/zz*/snpEffectPredictor.bin`
do
	DIR=`dirname $d`
	GEN=`basename $DIR`
	
	echo $GEN
	ZIP="snpEff_v"$SNPEFF_VERSION"_"$GEN".zip"
	zip -r $ZIP data/$GEN/*.bin
done

# Look for missing genomes
echo Missing genomes:
ls -d data/*/snpEffectPredictor.bin | grep -v genomes | cut -f 2 -d / | sort > genomes_bins.txt
ls -d data/* | grep -v genomes | cut -f 2 -d / | sort > genomes_dirs.txt
diff genomes_dirs.txt genomes_bins.txt | grep "^<"

