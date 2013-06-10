package ca.mcgill.mcb.pcingola.geneSets.algorithm;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import org.apfloat.Apfloat;

import ca.mcgill.mcb.pcingola.geneSets.GeneSet;
import ca.mcgill.mcb.pcingola.geneSets.GeneSets;
import ca.mcgill.mcb.pcingola.geneSets.Result;

/**
 * A generic enrichment algorithm for selecting gene-sets from a collection of gene-sets
 * 
 * @author pcingola
 */
public abstract class EnrichmentAlgorithm {

	public static final int HTML_TD_GENES_MAX_LEN = 40;
	public static final String HTML_BG_COLOR[] = { "dddddd", "eeeeee" };
	public static final String HTML_BG_COLOR_TITLE = "cccccc";

	boolean debug = false;
	boolean verbose = false;
	boolean htmlTable = false;
	int minGeneSetSize = 0;
	int maxGeneSetSize = Integer.MAX_VALUE;
	int numberToSelect;
	double maxPValue = 1.0;
	GeneSets geneSets;
	Set<String> filterOutputGeneSets;

	public static long PRINT_SOMETHING_TIME = 5000; // Print something every X milliseconds

	public EnrichmentAlgorithm(GeneSets geneSets, int numberToSelect) {
		this.geneSets = geneSets;
		this.numberToSelect = numberToSelect;
	}

	/**
	 * Showld we show this result or should the output be filtered?
	 * @param result
	 * @return
	 */
	protected boolean isShow(Result result) {
		// Filter by pValue
		if (result.getPvalueAdjusted() > maxPValue) return false;

		// Filter by Gene set name (show only if it is in this list)
		if (filterOutputGeneSets != null) {
			String gsName = result.getLatestGeneSet().getName();
			boolean show = filterOutputGeneSets.contains(gsName);
			if (verbose) System.err.println("\tFilter output list. Show geneSet '" + gsName + "' : " + show);
			return show;
		}

		// OK
		return true;
	}

	/**
	 * Print after all genesets are shown
	 */
	void printEnd() {
		if (htmlTable) System.out.println("</table>");
	}

	/**
	 * Show a result
	 * @param it : Iteration or Rank
	 * @param result
	 */
	void printResult(int it, Result result) {
		if ((result != null) && (result.getLatestGeneSet() != null)) {

			if (isShow(result)) {
				GeneSet geneSet = result.getLatestGeneSet();

				StringBuffer interestingGenes = new StringBuffer();
				ArrayList<String> genesList = new ArrayList<String>();
				genesList.addAll(geneSet.getInterestingGenes());
				Collections.sort(genesList);
				for (String gene : genesList)
					interestingGenes.append(gene + " ");

				if (htmlTable) {
					String descr = geneSet.getDescription();
					String name = geneSet.getName();

					if (descr.startsWith("http://")) {
						name = "<a href=\"" + descr + "\">" + name + "</a>";
						descr = "<a href=\"" + descr + "\">link</a>";
					}

					String intGenes = interestingGenes.toString();
					if (interestingGenes.length() > HTML_TD_GENES_MAX_LEN) intGenes = "<textarea rows=1 cols=" + HTML_TD_GENES_MAX_LEN + ">" + interestingGenes + "</textarea>";

					String bgcolor = HTML_BG_COLOR[it % 2];

					System.out.println("\t<tr bgcolor=\"" + bgcolor + "\"> <td nowrap>" + it //
							+ "</td>\t<td nowrap>" + String.format("%.2e", result.getPvalue().doubleValue()) //
							+ "</td>\t<td nowrap>" + String.format("%.2e", result.getPvalueAdjusted()) //
							+ "</td>\t<td nowrap>" + name //
							+ "</td>\t<td nowrap>" + geneSet.size() //
							+ "</td>\t<td nowrap>" + descr //
							+ "</td>\t<td nowrap>" + intGenes //
							+ "</td>\t<td nowrap>" + geneSet.rankSum() //
							+ "</td>\t</tr>");
				} else System.out.println("\t" + it //
						+ "\t" + result.getPvalue() //
						+ "\t" + result.getPvalueAdjusted() //
						+ "\t" + geneSet.getName() //
						+ "\t" + geneSet.size() //
						+ "\t" + geneSet.getDescription() //
						+ "\t" + result.getGeneSets() //
						+ "\t" + interestingGenes //
						+ "\t" + geneSet.rankSum() //
				);
			}
		} else System.out.println("\t" + it + "\tNULL");
	}

	/**
	 * Print before all genesets are shown
	 */
	void printTitle() {
		if (htmlTable) System.out.println("<table border=0> <tr bgcolor=\"" + HTML_BG_COLOR_TITLE + "\"> <th>Rank</th>\t<th>p-value</th>\t<th>p-value adj</th>\t<th>Latest result</th>\t<th>Size</th>\t<th>Description</th>\t<th>Interesting genes </th>\t<th> Score </th> </tr>");
		else if (verbose) System.out.println("\tIteration\tp-value\tp-value adj\tLatest result\tSize\tDescription\tResult\tInteresting genes");
	}

	/**
	 * Calculate the pValue for a given geneSet
	 * @param geneSetList
	 * @return
	 */
	abstract Apfloat pValue(GeneSet geneSet);

	/**
	 * Create a new gene set using all gene sets and calculate pValue
	 * @param geneSetList
	 * @return
	 */
	Apfloat pValue(List<GeneSet> geneSetList) {
		GeneSet newGeneSet = new GeneSet(geneSetList, geneSets);
		return pValue(newGeneSet);
	}

	/**
	 * Select the 'best' gene sets
	 * @return
	 */
	public Result select() {
		Result best = new Result();

		//---
		// Calculate pValues for each gene set mathing our criteria
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

	public void setFilterOutputGeneSets(Set<String> filterOutputGeneSets) {
		this.filterOutputGeneSets = filterOutputGeneSets;
	}

	public void setHtmlTable(boolean htmlTable) {
		this.htmlTable = htmlTable;
	}

	public void setMaxGeneSetSize(int maxGeneSetSize) {
		this.maxGeneSetSize = maxGeneSetSize;
	}

	public void setMaxPValue(double maxPValue) {
		this.maxPValue = maxPValue;
	}

	public void setMinGeneSetSize(int minGeneSetSize) {
		this.minGeneSetSize = minGeneSetSize;
	}

	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}

}
