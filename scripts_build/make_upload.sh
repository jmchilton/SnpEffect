#!/bin/sh 

rm -rvf snpEff_v3_2*.zip snpEff 
./scripts_build/make.sh 
./scripts_build/distro.sh 
mv snpEff_v3_2_core.zip snpEff_development.zip 
/batch/done.sh 
scp snpEff_development.zip pcingola,snpeff@frs.sourceforge.net:/home/frs/project/s/sn/snpeff/

