#-------------------------------------------------------------------------------
#
# Show simmulated values from Reactome + GTEx "circuit" simmulation
#
#																Pablo Cingolani
#-------------------------------------------------------------------------------

savePdf <- TRUE
maxPval <- 10E-30
minNameCount <- 50

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
rnames <- d$entityName

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

minColor = 1

if( savePdf ) { pdf(width=20, height=20); }

# Plot values
for( i in 1:length(rnames) ) {
	#if( length( grep("insulin", rnames[i], fixed=TRUE) ) > 0 ) {
	if(
		# (d$enityId[i] == 74695) || (d$enityId[i] == 165690) || (d$enityId[i] == 373676)
		#|| 
		( length( grep("insulin", rnames[i], fixed=TRUE) ) > 0 )
		) {
		x <- as.numeric(d[i,minExp:maxExp])

		# Kruskal-Wallis test
		kw <- kruskal.test( list( x, nameShort ) )
		pval.kw <- kw$p.value

		# Annova test
		anv <- oneway.test(x ~ nameShort)
		pval <- anv$p.value

		# Density plots
		if( !is.na(pval) && (pval < maxPval)) {
			name = as.character(d$entityName[i])

			# Boxplots
			par( mfrow=c(2,1) ) 
			boxplot(x ~ nameShort, main=name )

			cat('  \t', i,'\tp-value:', pval, '\tNode:', i, '\tId:', d$enityId[i], '\tName:', name, '\n');
			pvalStr <- paste('p-value(ANNOVA):', pval, '  p-value(Kruskal):', pval.kw)
			dens <- density(x)
			plot( dens, xlim=c(-1,1), main=rnames[i], sub=pvalStr )

			col <- minColor
			labels <- list()
			for( ns in nshort ) {
				keep <- (nameShort == ns) & (! is.na(x))

				xns <- x[keep]
				xns <- xns[ ! is.na(xns) ]
				if( length(xns) > minNameCount ) {
					lines( density(xns) , col=col )
					col <- col + 1
					labels[[ length(labels) + 1 ]] <- ns
					# cat('\t', ns, ':', length(xns), '\n');
				}
			}

			# Draw lables
			maxy <- max(dens$y)
			colors <- ( 1:length(labels) ) + minColor
			legend(0.75, maxy, labels, col=colors, lty=1 );
		} else {
			cat('NO\t', i,'\tp-value:', pval, '\tNode:', i, d$enityId[i], d$entityName[i],as.character(rnames[i]), '\n');
		}
	}
}

if( savePdf ) { dev.off(); }
