#!/bin/sh

# Get header from first file
head -n 1000 $1 | grep "^#"

# Concatenate all other files (without any header)
cat $* | grep -v "^#"

