#!/bin/sh

# ENSEMBL
#./scripts/queue.pl 20 20 15 queue_build.txt

# NCBI
./scripts/queue.pl 20 20 1 queue_build_ncbi.txt

grep "CDS check" *.stdout | cut -f 3- -d : | grep -v "^$" | sort | tee cds_check.txt
grep "Protein check" *.stdout | cut -f 3- -d : | grep -v "^$" | sort | tee protein_check.txt

 
