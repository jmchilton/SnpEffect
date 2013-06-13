package ca.mcgill.mcb.pcingola.geneSets.algorithm;

import java.util.HashMap;
import java.util.Random;

import org.apfloat.Apfloat;

import ca.mcgill.mcb.pcingola.geneSets.GeneSet;
import ca.mcgill.mcb.pcingola.geneSets.GeneSets;
import ca.mcgill.mcb.pcingola.gsa.PvalueList;
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
public class LeadingEdgeFractionAlgorithm extends EnrichmentAlgorithm {

	public static final double P_VALUE_CUTOFF_QUANTILE_DEFAULT = 0.95; // This it the value used in the paper
	public static final int FRACTIONS = 10000;
	public static final int MAX_FRACTIONS = 1000 * FRACTIONS;
	public static final Apfloat ONE = new Apfloat(1.0);

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
	PvalueList distribution(int size, double fraction) {
		PvalueList list = new PvalueList();

		int countFractions = 0;
		for (int i = 0; (i < FRACTIONS) && (countFractions > 0); i++) {
			double lef = randLeadingEdgeFraction(size);
			if (lef >= fraction) countFractions++;
			list.add(lef);

			if (i > MAX_FRACTIONS) break; // Too many iterations
		}

		return list;
	}

	/**
	 * Initialize GenePvalues
	 */
	void initGenePvalues() {
		GeneSet all = new GeneSet("", "", geneSets);
		for (String geneId : geneSets.getGenes())
			all.addGene(geneId);

		// Create an array of p-values
		genePvalues = new double[all.sizeEffective()];
		int i = 0;
		for (String geneId : geneSets.getGenes())
			if (geneSets.hasValue(geneId)) genePvalues[i++] = geneSets.getValue(geneId);
	}

	/**
	 * Calculate the pvalue of having a leading edge fraction if 'fraction' or more
	 * @param fraction
	 * @param size
	 * @return
	 */
	double leadingEdgeFractionPvalue(double fraction, int size) {
		if (size == 0) return 1.0;

		// Initialize genePvalues
		if (genePvalues == null) initGenePvalues();

		// Do we have a distribution?
		PvalueList pvals = distribution.get(size);
		if (pvals == null) {
			pvals = distribution(size, fraction);
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
		int count = 0, tot = 0;
		for (String geneId : geneSet) {
			if (geneSets.hasValue(geneId)) {
				if (geneSets.getValue(geneId) <= pValueCutOff) count++;
				tot++;
			}
		}

		// No genes have values? We are done
		if (tot <= 0) return ONE;

		// Calculate 'leading edge fraction'
		double leadingEdgeFraction = ((double) count) / ((double) tot);
		if (debug) Timer.showStdErr("Gene set: " + geneSet.getName() + "\tsize: " + geneSet.size() + "\tsize (eff): " + geneSet.sizeEffective() + "\t" + count + "\tleadingEdgeFraction: " + leadingEdgeFraction);

		// Calculate p-value
		double pvalue = leadingEdgeFractionPvalue(leadingEdgeFraction, tot);

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
		return lef;
	}

	public void setpValueCutOffQuantile(double pValueCutOffQuantile) {
		this.pValueCutOffQuantile = pValueCutOffQuantile;
	}
}
