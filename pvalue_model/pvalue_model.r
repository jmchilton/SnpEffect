
#-------------------------------------------------------------------------------
#
# Compare empirical distribution of ramdom reads to a Binomial 
# and Poisson estimates
#
#																Pablo Cingolani
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Histogram
#-------------------------------------------------------------------------------
histDens <- function( x, title, q=1.0, breaks = 50 ) {
    # Show only this part of the data
    xmin <- quantile( x, 1-q )
    xmax <- quantile( x, q )
    data <- x[ (x >= xmin) & (x <= xmax) ];

    dens <- density(data)

    h <- hist(data, main=title, xlab = "", ylab = "Frequency", freq = T, breaks=breaks);

    # Adjust density height to 'frecuency'
    dens$y <- max(h$counts) * dens$y/max(dens$y)
    lines(dens, col='red')

    # Mean & median calculated over the whola data
    #abline( v=mean(x), col='blue', lty=2, lwd=2);
    #abline( v=median(x), col='green', lty=2, lwd=2);

    #legend("topright",c("Mean","Median"),lty=c(1,1),col=c("blue","green"))
	return(h);
}

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

readSize <- 100
save <- TRUE

if( save )	png( width=1536, height=1024 );

# File having random samples (empirical distributions of random reads)
sampleFiles <- c( 'samples_10K_100_1K.txt', 'samples_10K_100_10K.txt',  'samples_10K_100_100K.txt' );  

# For each data file
for( sampleFile in sampleFiles ) {
	cat('\nFile:', sampleFile,'\n');

	# Load data
	samples <- read.table(sampleFile, sep="\t", header=T);
	markers <- read.table("markers.txt", sep="\t", header=T);

	reads <- samples$Chromosome[1]		# All reads must hit a chromosome, so this number is the number of reads

	# There is a rule of thumb stating that the Poisson distribution is a good approximation 
	# of the binomial distribution if n is at least 20 and p is smaller than or equal to 0.05, and 
	# an excellent approximation if n ≥ 100 and np ≤ 10
	chrMarker <- markers[ markers$type=='Chromosome',]	# Genome total size
	markers$p <- (markers$size + (readSize) * markers$count) / (chrMarker$size - readSize * chrMarker$count );	# Probability of a single Bernoulli trial
	markers$lambda <- markers$p * reads;		# Poisson's Lambda
	markers$poisson <- markers$lambda < 10		# Can we use Poisson's model?

	# Show histogram
	keep <- (markers$type != 'Chromosome') & (markers$type != 'Cds') & (markers$type != 'RareAminoAcid')
	par( mfrow=c(3,4) )

	cat('type\tmedian [real]\tmedian [model]\tstdDev [real]\tstdDev [model]\n')
	for( type in sort(markers$type[keep]) ) {
		thisMarker <- markers[ markers$type==type,]
		numReads <- samples[ ,names(samples) == type]

		# Show histogram (keep values in order to scale probabilities to frequencies)
		h <- histDens(numReads, type)

		# Show Poissin  model
		x <- min(numReads):max(numReads)
		p <- dpois(x, thisMarker$lambda)
		f <- max(h$counts) * p / max(p)	# Scale prob. to frequencies
		lines( x, f, col='blue')

		# Show binomial model
		p <- dbinom(x, reads, thisMarker$p)
		f <- max(h$counts) * p / max(p)	# Scale prob. to frequencies
		lines( x, f, col='green')

		# Show empirical vs estimated values
		binVar <- reads * thisMarker$p * (1 - thisMarker$p)
		cat(type, '\t', median(numReads), '\t', qbinom(0.5, reads, thisMarker$p), '\t', sqrt(var(numReads)), '\t', sqrt(binVar),'\n');
	}
}

if( save )	dev.off();
