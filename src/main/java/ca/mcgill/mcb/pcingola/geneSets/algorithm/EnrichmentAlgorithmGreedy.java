package ca.mcgill.mcb.pcingola.geneSets.algorithm;

import java.util.Date;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;

import org.apfloat.Apcomplex;
import org.apfloat.Apfloat;

import ca.mcgill.mcb.pcingola.geneSets.GeneSet;
import ca.mcgill.mcb.pcingola.geneSets.GeneSets;
import ca.mcgill.mcb.pcingola.geneSets.Result;

/**
 * A generic greedy enrichment algorithm for selecting gene-sets from a collection of gene-sets
 * 
 * @author pcingola
 */
public abstract class EnrichmentAlgorithmGreedy extends EnrichmentAlgorithm {

	public static final double DEFAULT_MAX_PVALUE = 0.05;

	protected boolean adjustedPvalue = true;
	protected double maxPvalue = DEFAULT_MAX_PVALUE;
	protected double maxPvalueAjusted = DEFAULT_MAX_PVALUE;

	public EnrichmentAlgorithmGreedy(GeneSets geneSets, int numberToSelect) {
		super(geneSets, numberToSelect);
	}

	/**
	 * Calculate best list of terms by adding a new term to a list that minimize p-value (rank sum)
	 * @return
	 */
	Result greedyPvalue(Result prevResult) {
		Apfloat pValue = Apcomplex.ONE;
		int geneSetCount = 0;
		Result best = new Result(prevResult.getList(), 1, 0);
		HashSet<GeneSet> genesetSet = new HashSet<GeneSet>();
		if (prevResult.getList() != null) genesetSet.addAll(prevResult.getList());
		Date start = new Date(), latest = new Date();

		// For each term...
		for (GeneSet geneSet : geneSets) {
			// Check GeneSet's conditions
			if ((geneSet.getGeneCount() > 0) // This term is empty? => skip it
					&& ((genesetSet == null) || (!genesetSet.contains(geneSet))) // Is this term already in the list? => skip it
					&& (geneSet.getGeneCount() >= minGeneSetSize) // Use gene sets bigger than minGeneSetSize
					&& (geneSet.getGeneCount() <= maxGeneSetSize) // Use gene sets smaller than maxGeneSetSize
			) {
				// Create a list of terms by joining the original list and adding a new term
				List<GeneSet> geneSetListNew = new LinkedList<GeneSet>();
				if (genesetSet != null) geneSetListNew.addAll(genesetSet);
				geneSetListNew.add(geneSet);

				// Calculate p-value
				pValue = pValue(geneSetListNew);

				// Is it better? => Store it
				if ((pValue.compareTo(Apcomplex.ZERO) > 0) && (pValue.compareTo(best.getPvalue()) < 0)) {
					// Copy latest list
					LinkedList<GeneSet> list = new LinkedList<GeneSet>();
					list.addAll(geneSetListNew);
					best.setPvalue(pValue);
					best.setList(list);
				}

				// Show something every now and then?
				Date now = new Date();
				long elapsed = now.getTime() - latest.getTime();
				long elapsedStart = now.getTime() - start.getTime();
				if (verbose && (elapsed > PRINT_SOMETHING_TIME)) {
					latest = now;
					System.err.println("\t\t\tElapsed:" + (elapsedStart / 1000) + " secs\tGene sets: " + geneSetListNew + "\tpValue: " + pValue + "\tbestPvalue: " + best.getPvalue() + "\t" + best.getList());
				}

				geneSetCount++;
			}
		}

		best.setGeneSetCount(geneSetCount); // This is used in order to adjust pValue
		return best;
	}

	@Override
	void printTitle() {
		if (htmlTable) System.out.println("<table border=0> <tr bgcolor=\"" + HTML_BG_COLOR_TITLE + "\"> <th>Iteration</th>\t<th>p-value</th>\t<th>p-value adj</th>\t<th>Latest result</th>\t<th>Size</th>\t<th>Description</th>\t<th>Interesting genes </th>\t<th> Score </th> </tr>");
		else if (verbose) System.out.println("\tIteration\tp-value\tp-value adj\tLatest result\tSize\tDescription\tResult\tInteresting genes");
	}

	/**
	 * Select the 'best' gene sets
	 * @return
	 */
	@Override
	public Result select() {
		printTitle();

		Result result = new Result();
		int iteration;
		for (iteration = 1; iteration <= numberToSelect; iteration++) {
			result = greedyPvalue(result);
			if (verbose) printResult(iteration, result); // Show something
			if (stopCriteria(result)) {
				if (verbose) System.out.println("\tStop criteria met.");
				break;
			}
		}

		if (htmlTable) System.out.println("</table>");
		else if (!verbose) printResult(iteration - 1, result);

		return result;
	}

	public void setAdjustedPvalue(boolean adjustedPvalue) {
		this.adjustedPvalue = adjustedPvalue;
	}

	@Override
	public void setMaxGeneSetSize(int maxGeneSetSize) {
		this.maxGeneSetSize = maxGeneSetSize;
	}

	public void setMaxPvalue(double maxPvalue) {
		this.maxPvalue = maxPvalue;
	}

	public void setMaxPvalueAjusted(double maxPvalueAjusted) {
		this.maxPvalueAjusted = maxPvalueAjusted;
	}

	@Override
	public void setMinGeneSetSize(int minGeneSetSize) {
		this.minGeneSetSize = minGeneSetSize;
	}

	@Override
	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}

	/**
	 * Stop criteria
	 * @param result
	 * @return true if stop criteria has been met and algorithm should stop iterating.
	 */
	protected boolean stopCriteria(Result result) {
		// No result? => Stop
		if (result == null) return true;

		// No geneSet selected? => Stop
		GeneSet geneSet = result.getLatestGeneSet();
		if (geneSet == null) return true;

		// Compare p-value to 'maxPvalue'
		if (adjustedPvalue) return result.getPvalueAdjusted() > maxPvalueAjusted;
		return result.getPvalue().doubleValue() > maxPvalue;
	}

}
