#!/bin/sh

cat zzz.500.txt | ../../../scripts_build/reactome/filterAllZero.pl > z.txt
cut -f 1,2 z.txt | tr -d "-" | tr -s "_" > ids.txt
cut -f 3- z.txt > nums.txt

wc -l nums.txt ids.txt z.txt

# Create file
head -n 1 zzz.500.txt |  tr -d "-" | tr -s "_" > zz.txt
paste ids.txt nums.txt >> zz.txt

head -n 1 zz.txt | tr "\t" "\n" > expNames.long.txt
head -n 1 zz.txt | tr "\t" "\n" | cut -f 1 -d "_" > expNames.short.txt
( echo -e "nameLong\tnameShort" ; paste expNames.long.txt expNames.short.txt ) > expNames.txt

#rm -vf ids.txt nums.txt
