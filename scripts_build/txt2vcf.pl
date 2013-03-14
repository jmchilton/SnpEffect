#!/usr/bin/perl

#-------------------------------------------------------------------------------
#
# Convert simple TXT files into VCF-like files
#
#
# 															Pablo Cingolani
#-------------------------------------------------------------------------------

for( $lineNum=1 ; $l = <STDIN> ; $lineNum++ ) {
	chomp $l;
	$l =~ tr/\n\r//d;

	($chr, $pos, $ref, $alt, $val) = split /\t/, $l;

	# Header? Use value as INFO field name
	if( $lineNum == 1 ) {
		$infoName = $val;		
		print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
	} else {
		print "$chr\t$pos\t.\t$ref\t$alt\t.\t.\t$infoName=$val\n";
	}
	
}
