#!/bin/sh

SNPEFF="java -Xmx2G -jar snpEff.jar"

# Test cases hg37
$SNPEFF build -v -txt testCase

# Test cases hg37.61
$SNPEFF build -v -gtf22 testHg3761Chr15
$SNPEFF build -v -gtf22 testHg3761Chr16 

# Test cases hg37.63
$SNPEFF build -v -gtf22 testHg3763Chr1
$SNPEFF build -v -gtf22 testHg3763Chr20 
$SNPEFF build -v -gtf22 testHg3763ChrY 

# Test cases hg37.65
$SNPEFF build -v -gtf22 testHg3765Chr22
