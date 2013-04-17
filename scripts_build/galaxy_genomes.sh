#!/bin/sh

# Databases list
# Remove the first two lines (titles)
(
	echo "# SnpEff databases"
	echo "# List created using command: java -jar snpEff.jar databases"
	echo "#Version    Description"
	java -jar snpEff.jar databases | tail +3 | cut -f 1,2 
) > galaxy/tool-data/snpeffect_genomedb.loc.sample
