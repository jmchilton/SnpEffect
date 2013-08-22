#-------------------------------------------------------------------------------
#
# P-values correction for gene set analysis (SnpEff gsa)
#
#
#															Pablo Cingolani 2013
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

show = FALSE
# show = TRUE		# Used for debugging

#---
# Parse command line arguments
#---
if( ! exists('fileName') ) {
    args <- commandArgs(trailingOnly = TRUE)
    fileName = args[1]
    outFileName = args[2]
}

#---
# Read data
#---
d = read.table(fileName, sep="\t", header=TRUE)

p = d$score			# p-values or scores
c = d$scoreCount	# number of variants 

keep = (p > 0) & (c > 0)	# Filter out p-values of zero OR counts of zero (because of log)
p = p[keep]
lp = qnorm(p)				# Convert to z-scores using an inverse normal distribution

c = c[keep]
lc = log10(c)

#---
# Linear model
#---
lmfit = lm( lp ~ lc )
sumLmfit =  summary(lmfit)
pvalCoef = sumLmfit$coefficients[2,4]

res = lmfit$residuals		# Residuals
padj = pnorm(res)			# Adjusted p-values

if( show ) {
	print(sumLmfit)
	cat('Slope:\t', sumLmfit$coefficients[2], '\tp-value:\t', pvalCoef, 'File:\t', fileName, '\n');

	par(mfrow=c(2,2))
	smoothScatter( lc, lp, main="Z-Scores", xlab="Number of scores", ylab="Z-Score" )
	lines(lowess(lc, lp ), col='orange')

	smoothScatter( lc, res, main="Adjusted Z-Scores", xlab="Number of scores", ylab="Z-Score" )
	lines(lowess(lc, res ), col='orange')

	plot( density(lp) , main="Z-Scores distribution", xlab="red:Adjusted black:Unadjusted")
	lines( density(res) , col='red')

	plot( density(lc) , main="log10(Number of scores)", xlab="")
}

#---
# Decide whether to use corrected scores or not
#---
if( pvalCoef < 10^-12 ) { 
	so = padj		# Significant? Then use correction
} else {
	so = d$score	# Not very significant? Use original scores
}

#---
# Create output file
#---
dout = data.frame( geneId = d$geneId[keep], score = so )
write.table( dout, file=outFileName, quote=FALSE, row.names = FALSE, col.names = FALSE, sep="\t")

