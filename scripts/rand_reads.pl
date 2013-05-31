#!/usr/bin/perl


$readSize = 50;
$chrSize = 23011544;
$numReads = 1 * 1000 * 1000;

for( $i=0 ; $i < $numReads ; $i++ ) {
	$start = int( rand() * $chrSize );
	$end = $start + $readSize;
	print "2L\t$start\t$end\n";
}
