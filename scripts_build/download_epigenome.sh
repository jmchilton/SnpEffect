#!/bin/sh -e

# Download current data release from Epigenome project

#wget -r --no-parent http://www.genboree.org/EdaccData/Current-Release/sample-experiment/Pancreatic_Islets/
#wget -r --no-parent http://www.genboree.org/EdaccData/Current-Release/sample-experiment/Adult_Liver/
#wget -r --no-parent http://www.genboree.org/EdaccData/Current-Release/sample-experiment/Adipose_Tissue/
#wget -r --no-parent http://www.genboree.org/EdaccData/Current-Release/sample-experiment/Skeletal_Muscle/

wget -nv -r --no-parent -A "*.bed.gz" http://www.genboree.org/EdaccData/Current-Release/sample-experiment/Skeletal_Muscle/
mv www.genboree.org/EdaccData/Current-Release/sample-experiment

