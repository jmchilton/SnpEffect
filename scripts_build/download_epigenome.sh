#!/bin/sh

cd db/epigenome/

#wget -r --no-parent --no-clobber -A "*.bed.gz" http://www.genboree.org/EdaccData/Current-Release/sample-experiment/Pancreatic_Islets/
#wget -r --no-parent --no-clobber -A "*.bed.gz" http://www.genboree.org/EdaccData/Current-Release/sample-experiment/Adult_Liver/
#wget -r --no-parent --no-clobber -A "*.bed.gz" http://www.genboree.org/EdaccData/Current-Release/sample-experiment/Adipose_Tissue/
#wget -r --no-parent --no-clobber -A "*.bed.gz" http://www.genboree.org/EdaccData/Current-Release/sample-experiment/Skeletal_Muscle/

url="http://www.genboree.org/EdaccData/Current-Release/sample-experiment"

# Download all antries in directory
for dir in `wget -q -O - $url | cut -f 2 -d \" | grep "/$" | grep -v "\.\./"`
do
    echo $dir
    wget -nv -r --no-parent --no-clobber -A "*.bed.gz" "$url/$dir"
done

