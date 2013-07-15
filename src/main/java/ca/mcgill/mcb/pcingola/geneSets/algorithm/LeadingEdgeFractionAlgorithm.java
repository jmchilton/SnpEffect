package ca.mcgill.mcb.pcingola.geneSets.algorithm;

import org.apfloat.Apfloat;

import ca.mcgill.mcb.pcingola.geneSets.GeneSet;
import ca.mcgill.mcb.pcingola.geneSets.GeneSets;
import ca.mcgill.mcb.pcingola.gsa.ScoreList;
import ca.mcgill.mcb.pcingola.probablility.FisherExactTest;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * Leading edge fraction algorithm
 * 
 * References: "Common Inherited Variation in Mitochondrial Genes Is Not Enriched for Associations with Type 2 Diabetes or Related Glycemic Traits"
 * 				http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1001058 
 * 				See page 12, "Step 4"
 * 
 * @author pablocingolani
 */
public class LeadingEdgeFractionAlgorithm extends FisherPValueGreedyAlgorithm {

	public static final double P_VALUE_CUTOFF_QUANTILE_DEFAULT = 0.95; // This it the value used in the paper

	public static final Apfloat ONE = new Apfloat(1.0);

	double pValueCutOff;
	double pValueCutOffQuantile = P_VALUE_CUTOFF_QUANTILE_DEFAULT;
	int N, D;

	public LeadingEdgeFractionAlgorithm(GeneSets geneSets, int numberToSelect) {
		super(geneSets, numberToSelect);
		init();
	}

	/**
	 * Initialize parameters
	 */
	void init() {
		// Calculate pvalue cutoff
		pValueCutOff = pValueCutOff();

		//---
		// Calculate Fisher parameters
		//---
		N = D = 0;
		for (String geneId : geneSets.getGenes())
			if (geneSets.hasValue(geneId)) {
				N++;
				if (geneSets.getValue(geneId) <= pValueCutOff) D++;
			}

		if (verbose) Timer.showStdErr("Fisher Exact test parameters:\n\tN : " + N + "\n\tD : " + D);
	}

	/**
	 * Create a new gene set using all gene sets and calculate pValue
	 * @param geneSetList
	 * @return
	 */
	@Override
	Apfloat pValue(GeneSet geneSet) {
		// Count number of p-values less or equal to 'pValueCutOff'
		int count = 0, tot = 0;
		for (String geneId : geneSet) {
			if (geneSets.hasValue(geneId)) {
				if (geneSets.getValue(geneId) <= pValueCutOff) count++;
				tot++;
			}
		}

		// No genes have values? We are done
		if (tot <= 0) return ONE;

		if (debug) {
			// Calculate and show 'leading edge fraction'
			double leadingEdgeFraction = ((double) count) / ((double) tot);
			Timer.showStdErr("Gene set: " + geneSet.getName() + "\tsize: " + geneSet.size() + "\tsize (eff): " + geneSet.sizeEffective() + "\t" + count + "\tleadingEdgeFraction: " + leadingEdgeFraction);
		}

		// Calculate p-value
		double pvalueFisher = FisherExactTest.get().fisherExactTestUp(count, N, D, tot);
		return new Apfloat(pvalueFisher);
	}

	/**
	 * Calculate 'pValueCutOff' (see paper methods)
	 * @return
	 */
	double pValueCutOff() {
		// Create a list of p-values
		ScoreList pvlist = new ScoreList();
		for (String geneId : geneSets.getGenes())
			if (geneSets.hasValue(geneId)) pvlist.add(geneSets.getValue(geneId));

		double pco = pvlist.quantile(1 - pValueCutOffQuantile);

		// Show
		if (verbose) Timer.showStdErr("Calculate pValue_CutOff: " //
				+ "\n\tSize (effective) : " + pvlist.size() //
				+ "\n\tQuantile         : " + pValueCutOffQuantile //
				+ "\n\tp-value CutOff   : " + pco //
		);
		if (debug) Timer.showStdErr("\tp-values: " + pvlist);

		return pco;
	}

	public void setpValueCutOffQuantile(double pValueCutOffQuantile) {
		this.pValueCutOffQuantile = pValueCutOffQuantile;
		pValueCutOff = pValueCutOff();
	}
}
