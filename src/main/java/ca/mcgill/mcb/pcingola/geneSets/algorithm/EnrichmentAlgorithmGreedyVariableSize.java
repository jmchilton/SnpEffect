package ca.mcgill.mcb.pcingola.geneSets.algorithm;

import ca.mcgill.mcb.pcingola.geneSets.GeneSets;

/**
 * A greedy enrichment algorithm for selecting gene-sets using a variable geneSet-size strategy:
 * 
 * 	i) Select only from geneSets in low-sizes e.g. geneSet.size() in [1-10]
 * 	ii) If p-value goes up, use a larger geneSize range, e.g. geneSet.size() in [11-20]
 * 
 * @author pcingola
 */
public abstract class EnrichmentAlgorithmGreedyVariableSize extends EnrichmentAlgorithmGreedy {

	public EnrichmentAlgorithmGreedyVariableSize(GeneSets geneSets, int numberToSelect) {
		super(geneSets, numberToSelect);
	}

}
