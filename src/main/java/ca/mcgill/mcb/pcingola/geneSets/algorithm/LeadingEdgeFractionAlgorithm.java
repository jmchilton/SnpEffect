package ca.mcgill.mcb.pcingola.geneSets.algorithm;

import java.util.HashMap;
import java.util.Random;

import org.apfloat.Apfloat;

import ca.mcgill.mcb.pcingola.geneSets.GeneSet;
import ca.mcgill.mcb.pcingola.geneSets.GeneSets;
import ca.mcgill.mcb.pcingola.gsa.PvalueList;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * Leading edge fraction algorithm
 * 
 * References: "Common Inherited Variation in Mitochondrial Genes Is Not Enriched for Associations with Type 2 Diabetes or Related Glycemic Traits"
 * 				http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1001058 
 * 				See page 12, "Step 4"
 * 
 * @author pablocingolani
 */
public class LeadingEdgeFractionAlgorithm extends EnrichmentAlgorithm {

	public static final double P_VALUE_CUTOFF_QUANTILE_DEFAULT = 0.95; // This it the value used in the paper
	public static final int FRACTIONS = 10000;

	Random random = new Random();
	double pValueCutOff = Double.NaN;
	double pValueCutOffQuantile = P_VALUE_CUTOFF_QUANTILE_DEFAULT;
	double genePvalues[];
	HashMap<Integer, PvalueList> distribution;

	public LeadingEdgeFractionAlgorithm(GeneSets geneSets, int numberToSelect) {
		super(geneSets, numberToSelect);
		distribution = new HashMap<Integer, PvalueList>();
	}

	/**
	 * Calculate a distribution of leading edge fractions of size 'size'
	 * @return
	 */
	PvalueList distribution(int size) {
		Gpr.debug("Calculating leading edge distribution for size : " + size);
		PvalueList list = new PvalueList();

		for (int i = 0; i < FRACTIONS; i++) {
			double lef = randLeadingEdgeFraction(size);
			list.add(lef);
		}

		return list;
	}

	/**
	 * Calculate the pvalue of having a leading edge fraction if 'fraction' or more
	 * @param fraction
	 * @param size
	 * @return
	 */
	double leadingEdgeFractionPvalue(double fraction, int size) {
		// Initialize genePvalues
		if (genePvalues == null) {
			//!!!!!!! WHAT HAPPENS WHEN MOST GENES DO NOT HAVE VALUES?
			genePvalues = new double[geneSets.getGenes().size()];
			// Populate with p-values
			int i = 0;
			for (String geneId : geneSets.getGenes())
				genePvalues[i++] = geneSets.getValue(geneId);
		}

		// Do we have a distribution?
		PvalueList pvals = distribution.get(size);
		if (pvals == null) {
			pvals = distribution(size);
			distribution.put(size, pvals);
		}

		return 1.0 - pvals.cdf(fraction);
	}

	/**
	 * Create a new gene set using all gene sets and calculate pValue
	 * @param geneSetList
	 * @return
	 */
	@Override
	Apfloat pValue(GeneSet geneSet) {
		if (Double.isNaN(pValueCutOff)) pValueCutOff = pValueCutOff();

		// Count number of pvalues less or equal to 'pValueCutOff'
		int count = 0;
		for (String geneId : geneSet)
			if (geneSets.getValue(geneId) <= pValueCutOff) count++;

		// Claculate 'leading edge fraction'
		double leadingEdgeFraction = ((double) count / ((double) geneSet.size()));

		// Calculate p-value
		double pvalue = leadingEdgeFractionPvalue(leadingEdgeFraction, geneSet.size());

		return new Apfloat(pvalue);
	}

	/**
	 * Calculate 'pValueCutOff' (see paper methods)
	 * @return
	 */
	double pValueCutOff() {
		// Create a list of p-values
		PvalueList pvlist = new PvalueList();
		for (String geneId : geneSets.getGenes())
			pvlist.add(geneSets.getValue(geneId));

		return pvlist.quantile(1 - pValueCutOffQuantile);
	}

	/**
	 * Randomly sample 'size' p-values and calculate leading edge fraction 
	 * @param size
	 * @return
	 */
	double randLeadingEdgeFraction(int size) {
		int count = 0;
		for (int i = 0; i < size; i++) {
			int idx = random.nextInt(genePvalues.length);
			if (genePvalues[idx] <= pValueCutOff) count++;
		}

		double lef = ((double) count / size);
		Gpr.debug("\t" + lef);
		return lef;
	}

	public void setpValueCutOffQuantile(double pValueCutOffQuantile) {
		this.pValueCutOffQuantile = pValueCutOffQuantile;
	}
}
