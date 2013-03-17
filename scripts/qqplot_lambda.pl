#!/usr/bin/perl
#-------------------------------------------------------------------------------
#
# Plot a QQ plot (using R) 
# Data is feed as a 1 column of numbers 
#
# Note: Any line that does not match a numeric regular expression, is filtered out).
#
#														Pablo Cingolani
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

# Parse command line option (file base name)
$base = 'QQ-plot';
$plot = 1;
if( $ARGV[0] ne '' )	{ $base = $ARGV[0]; }
if( $ARGV[1] eq '-noPlot' )	{ $plot = 0; }

$pngFile = "$base.png";
$txtFile = "$base.txt";

# Read STDIN and create an R vector
open TXT, "> $txtFile" or die "Cannot open output file '$txtFile'\n";
print TXT "x\n";
for( $ln = 0 ; $l = <STDIN> ; ) {
	chomp $l;

	# Does the string contain exactly one number? (can be float)
	if( $l =~ /^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/ ) { print TXT "$l\n"; }
}
close TXT;

#---
# Create an R program, save QQ-plot as PNG image
#---

open R, "| R --vanilla --slave " or die "Cannot open R program\n";
print R <<EOF;

library('GenABEL')

qqplot <- function( x, title ) {
	keep <- (x > 0) & (x <= 1) & ( ! is.na(x) );
	x <- x[keep]
	l <- estlambda(x, plot = F)
	cat('\nTitle:\t', title, '\tLambda:\t', l\$estimate, '\n');

	title <- paste( title , " lambda=", l\$estimate)
	estlambda(x, plot = TRUE, main=title)
}

png('$pngFile', width = 1024, height = 1024);

data <- read.csv("$txtFile", sep='\t', header = TRUE);
qqplot( data\$x, "$base" );

dev.off();
quit( save='no' )
EOF

close R;

#---
# Show figure
#---

if( $plot ) {
	$os = `uname`;
	$show = "eog"; 
	if( $os =~ "Darwin" )	{ $show = "open"; }
	`$show $pngFile &`;
}

