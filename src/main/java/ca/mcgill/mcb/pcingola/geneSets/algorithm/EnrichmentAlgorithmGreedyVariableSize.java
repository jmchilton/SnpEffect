package ca.mcgill.mcb.pcingola.geneSets.algorithm;

import ca.mcgill.mcb.pcingola.geneSets.GeneSet;
import ca.mcgill.mcb.pcingola.geneSets.GeneSets;
import ca.mcgill.mcb.pcingola.geneSets.Result;

/**
 * A greedy enrichment algorithm for selecting gene-sets using a variable geneSet-size strategy:
 * 
 * 	i) Select only from geneSets in low-sizes e.g. geneSet.size() in [1-10]
 * 	ii) If p-value goes up, use a larger geneSize range, e.g. geneSet.size() in [11-20]
 * 
 * @author pcingola
 */
public abstract class EnrichmentAlgorithmGreedyVariableSize extends EnrichmentAlgorithmGreedy {

	int sizeGrowthFactor = 2;
	int initialSize = 50;
	int maxSize = -1;

	public EnrichmentAlgorithmGreedyVariableSize(GeneSets geneSets, int numberToSelect) {
		super(geneSets, numberToSelect);

		// Update maxSize
		for (GeneSet gs : geneSets)
			maxSize = Math.max(maxSize, gs.size());
	}

	/**
	 * Find best p-value using a greedy algorithm (smaller gene sets first)
	 */
	@Override
	protected Result greedyPvalue(Result prevResult) {
		minGeneSetSize = 1;
		maxGeneSetSize = initialSize;

		// Focus on smaller gene set sizes first
		while (minGeneSetSize < maxSize) {
			if (debug) System.err.println("Trying size range [ " + minGeneSetSize + " , " + maxGeneSetSize + " ]");

			// Run greedy algorithm from parent class
			Result best = super.greedyPvalue(prevResult);

			// Update real count size.
			// Note: We are changing the search order, but "search space size" remains the same as when we have all geneSets.
			best.setGeneSetCountLast(geneSets.getGeneSetCount() - best.getGeneSets().size());

			// Found something better? => we are done
			double bestPval = best.getPvalueAdjusted();
			if (bestPval <= prevResult.getPvalueAdjusted() && (bestPval <= maxPvalueAjusted)) return best;

			// Update new size range
			minGeneSetSize = maxGeneSetSize + 1;
			maxGeneSetSize = sizeGrowthFactor * maxGeneSetSize;
		}

		return null;
	}

	public void setInitialSize(int initialSize) {
		this.initialSize = initialSize;
	}

	public void setSizeGrowthFactor(int sizeGrowthFactor) {
		this.sizeGrowthFactor = sizeGrowthFactor;
	}
}
