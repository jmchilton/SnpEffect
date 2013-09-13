
maxPval <- 10E-30
minNameCount <- 10

if( ! exists('d') ) {
	# Simaltion data
	d <- read.table('zz.txt', sep='\t', header=TRUE)

	# Experiment names
	names <- read.table('expNames.txt', sep='\t', header=TRUE)
}

# Experiment indexes
minExp <- 3
maxExp <- dim(d)[2]

# Column names
cnames <- colnames(d)
rnames <- rownames(d)

# Experiment name
# Remove low counts and 'NA'
nameShort <- as.vector( names$nameShort[minExp:maxExp] )
nameShort[ is.na(nameShort) ] <- "Etc"
nshort <- sort(unique(nameShort))
for( ns in nshort ) {
	keep <- (nameShort == ns)
	if( sum(nameShort == ns) <= minNameCount ) { 
		nameShort[ nameShort == ns ] <- "Etc" 
	}
}

# Plot values
for( i in 1:length(rnames) ) {
	x <- as.numeric(d[i,minExp:maxExp])


	# Kruskal-Wallis test
	kw <- kruskal.test( list( x, nameShort ) )
	pval.kw <- kw$p.value

	# Annova test
	anv <- oneway.test(x ~ nameShort)
	pval <- anv$p.value

	if( is.na(pval) && !is.na(pval.kw) ) {
		pval <- pval.kw
	} else if( !is.na(pval) && !is.na(pval.kw) && (pval < pval.kw)) { 
		cat("KW!\n") 
		pval <- max(pval, pval.kw);
	}

	if( !is.na(pval) && (pval < maxPval)) {
		cat('  \t', i,'\tp-value:', pval, '\tNode:', i, rnames[i], '\n');
		plot( density(x), main=rnames[i] )

		col <- 0
		for( ns in nshort ) {
			keep <- (names$nameShort == ns) & (! is.na(x))

			xns <- x[keep]
			xns <- xns[ ! is.na(xns) ]
			if( length(xns) > minNameCount ) {
				lines( density(xns) , col=col )
			}
			col <- col + 1
		}
	} else {
		cat('NO\t', i,'\tp-value:', pval, '\tNode:', i, rnames[i], '\n');
	}
}
