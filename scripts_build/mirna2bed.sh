#!/bin/sh

cat db/miRNA/human_predictions*.txt \
	| ./scripts_build/mirna2bed.pl \
	| sort -k1,1 -k2,2g \
	> db/miRNA/human_predictions.bed

bgzip human_predictions.bed
tabix -p bed human_predictions.bed.gz
