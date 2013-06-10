package ca.mcgill.mcb.pcingola.geneSets.algorithm;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.apfloat.Apfloat;

import ca.mcgill.mcb.pcingola.geneSets.GeneSet;
import ca.mcgill.mcb.pcingola.geneSets.GeneSets;
import ca.mcgill.mcb.pcingola.geneSets.Result;

/**
 * A greedy enrichment algorithm for selecting gene-sets 
 * 
 * @author pcingola
 */
public abstract class EnrichmentAlgorithmGreedyVariableSize extends EnrichmentAlgorithmGreedy {

	public EnrichmentAlgorithmGreedyVariableSize(GeneSets geneSets, int numberToSelect) {
		super(geneSets, numberToSelect);
	}

	/**
	 * Select the 'best' gene sets
	 * @return
	 */
	@Override
	public Result select() {
		Result best = new Result(null, 1, 0);

		//---
		// Calculate pValues
		//---
		List<Result> results = new ArrayList<Result>();
		for (GeneSet geneSet : geneSets) {
			if ((geneSet.getGeneCount() > 0) // This term is empty? => skip it
					&& (geneSet.getGeneCount() >= minGeneSetSize) // Use gene sets bigger than minGeneSetSize
					&& (geneSet.getGeneCount() <= maxGeneSetSize) // Use gene sets smaller than maxGeneSetSize
			) {
				// Calculate pValue
				Apfloat pValue = pValue(geneSet);
				Result result = new Result(geneSet, pValue, 0); // We'll update the geneSetCount later
				results.add(result);
			}
		}

		// Update geneSetCounts
		for (Result r : results)
			r.setGeneSetCount(results.size());

		//---
		// Show results
		//---
		if (htmlTable || verbose) {
			printTitle();

			// Sort by pValue
			Collections.sort(results);

			// Show them
			int rank = 1;
			for (Result r : results)
				printResult(rank++, r);

			printEnd();
		}

		return best;
	}

}
