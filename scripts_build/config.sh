#!/bin/sh

export VERSION_SNPEFF=3.2
export SUBVERSION_SNPEFF="a"
export VERSION_SNPSIFT=1.9
export VERSION_SNPSQL=0.2

export ENSEMBL_RELEASE=70
export ENSEMBL_BFMPP_RELEASE=18

# Version values using underscores ('3_2' instead of '3.2')
export SNPEFF_VERSION=`echo $VERSION_SNPEFF | tr "." "_"`
export SNPSIFT_VERSION=`echo $VERSION_SNPSIFT | tr "." "_"`
export SNPSQL_VERSION=`echo $VERSION_SNPSQL | tr "." "_"`

export SNPEFF_VERSION_REV=$SNPEFF_VERSION"$SUBVERSION_SNPEFF"
