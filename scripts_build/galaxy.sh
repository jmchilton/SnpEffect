#!/bin/sh

# Create config file for Galaxy
# Note: Replaces 'data_dir' to a relative path
cat snpEff.config | sed "s/^data_dir.*/data_dir = data/" > galaxy/snpEff.config

# Create genomes databases list
# Remove the first two lines (titles)
mkdir -p galaxy/tool-data
(
	echo "# SnpEff databases"
	echo "# List created using command: java -jar snpEff.jar databases"
	echo "#Version    Description"
	java -jar snpEff.jar databases | tail -n +3 | cut -f 1,2 
) > galaxy/tool-data/snpeffect_genomedb.loc
cp galaxy/tool-data/snpeffect_genomedb.loc galaxy/tool-data/snpeffect_genomedb.loc.sample
