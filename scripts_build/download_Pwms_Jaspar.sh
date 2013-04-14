#!/bin/sh

mkdir -p db/jaspar/
cd db/jaspar/

wget "http://jaspar.genereg.net/html/DOWNLOAD/all_data/matrix_only/matrix_only.txt"
gzip matrix_only.txt
mv matrix_only.txt pwms.bin

echo "File pwms.bin created"
