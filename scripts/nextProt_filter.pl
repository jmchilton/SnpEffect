#!/usr/bin/perl

#-------------------------------------------------------------------------------
#
# Find variants having at least one HIGH or MODERATE impact 
# variant AND one NEXT_PROT (HIGH) impact
#
#
#																Pablo Cingolani
#-------------------------------------------------------------------------------

# Debug mode?
$debug=0;

# Read STDIN
while( $l = <STDIN> ) {

	# Always show VCF headers
	if( $l =~ /^#/ ) {
		print $l;
	} else {
		# Split VCF entry
		@t = split /\t/, $l;

		# Get INFO field
		$info = $t[7];
		print "INFO: $info\n" if $debug;

		# Parse INFO field: Get EFF tag
		$eff = '';
		@infos = split /;/, $info;
		foreach $info ( @infos ) {
			if( $info =~ /^EFF=(.*)/ ) { $eff = $1; }
		}

		# Parse EFF tag
		if( $eff ne '' ) {
			@effs = split /,/, $eff;

			# Try to find one EFF that is at least MODERATE and one NEXT_PROT that is HIGH
			$effBest = $nextprotBest = '';
			foreach $f ( @effs ) {
				if(( $f !~/^NEXT_PROT/ ) && (( $f =~ /HIGH/ ) || ( $f =~ /MODERATE/ )) ) {
					print "\tEFF      : $f\n" if $debug;
					$effBest=$f;
				} elsif(( $f =~/^NEXT_PROT/ ) && ( $f =~ /HIGH/ )) {
					print "\tNEXT_PROT: $f\n" if $debug;
					$nextprotBest = $f;
				}
			}

			# Do we have at least one Moderate EFFect AND one High NEXT_PROT?
			if(($effBest ne '') && ($nextprotBest ne '')) { 
				print "$effBest\t$nextprotBest\t" if $debug; 
				print "$l"; 
				} 
		}
	}

}
