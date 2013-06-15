package ca.mcgill.mcb.pcingola.gsa;

import flanagan.analysis.Stat;
import gnu.trove.list.array.TDoubleArrayList;

/**
 * A list of pvalues (i.e. values in the arnge [0, 1])
 * 
 * @author pcingola
 */
public class PvalueList {

	public enum PvalueSummary {
		MIN, AVG, AVG10, FISHER_CHI_SQUARE, Z_SCORES, SIMES, BONFERRONI, FDR
	}

	public static final double SIGNIFICANCE_LEVEL_95 = 0.05;

	String geneId;
	TDoubleArrayList pValues;
	boolean sorted = false;

	/**
	 * Upper tail 1 - ChiSquareCDF(p)
	 * @param chiSquare
	 * @param nu
	 * @return
	 */
	public static double chiSquareCDFComplementary(double chiSquare, int nu) {
		if (nu <= 0) throw new IllegalArgumentException("The degrees of freedom [nu], " + nu + ", must be greater than zero");
		return Stat.incompleteGammaComplementary(nu / 2.0D, chiSquare / 2.0D);
	}

	public PvalueList() {
		pValues = new TDoubleArrayList();
	}

	/**
	 * Add a p-value to the list
	 * @param pvalue
	 */
	public void add(double pvalue) {
		if ((pvalue < 0) || (pvalue > 1)) throw new RuntimeException("p-value out of range: " + pvalue);
		pValues.add(pvalue);
	}

	/**
	 * Cumulative distribution function of p-values: 
	 * 
	 * 			P[ pValues <= p ]		(i.e. lower tail).
	 * 
	 * @param p
	 * @return
	 */
	public double cdf(double p) {
		if (size() <= 0) return 1;

		sort();
		int idx = pValues.binarySearch(p);
		if (idx < 0) idx = -(idx + 1); // If 'p' is not found, idx is (-insertion_point - 1);

		return ((double) idx) / size();
	}

	/**
	 * Cumulative distribution function of p-values: 
	 * 
	 * 			P[ pValues > p ]		(i.e. upper tail).
	 * 
	 * @param p
	 * @return
	 */
	public double cdfUpper(double p) {
		if (size() <= 0) return 1;

		sort();
		int idx = pValues.binarySearch(p);
		if (idx < 0) idx = -(idx + 1); // If 'p' is not found, idx is (-insertion_point - 1);

		return ((double) (size() - idx)) / size();
	}

	public String getGeneId() {
		return geneId;
	}

	public double getPvalue(int index) {
		return pValues.get(index);
	}

	/**
	 * Create a single pValue representing the gene
	 * @return
	 */
	public double pValue(PvalueSummary pvalueSummary) {
		switch (pvalueSummary) {
		case MIN:
			return pValueMin();

		case AVG:
			return pValueAvg();

		case AVG10:
			return pValueAvgTop(10);

		case FISHER_CHI_SQUARE:
			return pValueFisherChi2();

		case Z_SCORES:
			return pValueZScore();

		case SIMES:
			return pValueSimes();

		case BONFERRONI:
			return pValueBonferroni();

		case FDR:
			return pValueFdr(SIGNIFICANCE_LEVEL_95);

		default:
			throw new RuntimeException("Unimplemented method for summary '" + pvalueSummary + "'");
		}
	}

	/**
	 * Get average pvalue
	 * @return
	 */
	public double pValueAvg() {
		if (size() <= 0) return 1.0;

		double sum = 0.0;
		for (int i = 0; i < size(); i++)
			sum += getPvalue(i);

		return sum / size();
	}

	/**
	 * Get average pvalue
	 * @return
	 */
	public double pValueAvgTop(int topN) {
		if (size() <= 0) return 1.0;

		sort(); // Sort collection

		// Sum smallest values
		double sum = 0.0;
		int max = Math.min(size(), topN);
		for (int i = 0; i < max; i++)
			sum += getPvalue(i);

		return sum / max;
	}

	/**
	 * Minimum p-value corrected using Bonferroni
	 * @return
	 */
	public double pValueBonferroni() {
		if (pValues.size() <= 0) return 1.0;
		return Math.min(1.0, pValueMin() * pValues.size());
	}

	/**
	 * Combine p-values using FDR procedure
	 * 
	 * References: http://en.wikipedia.org/wiki/False_discovery_rate
	 * 
	 * @return A combined p-value
	 */
	public double pValueFdr(double alpha) {
		if (size() <= 0) return 1.0;

		// Count non-zero p-values (we treat zero p-values as errors, so we skip them) 
		int total = 0;
		for (int i = 0; i < size(); i++)
			if (getPvalue(i) > 0) total++;
		double tot = total;

		// No p-value? => We are done
		if (total <= 0) return 1.0;

		sort(); // Sort collection

		// Perform Simes's method
		int count = 0;
		double pFdrMax = 0.0;
		for (int i = 0; i < size(); i++) {
			double pvalue = getPvalue(i);

			// We treat zero p-values as errors, so we skip them 
			if (pvalue > 0) {
				count++;
				double pFdr = tot * pvalue / count;
				if ((pFdr <= alpha) || (pFdr < pFdrMax)) pFdrMax = pFdr;
			}
		}

		return pFdrMax;
	}

	/**
	 * Combine p-values using Fisher's method
	 * 
	 * References: http://en.wikipedia.org/wiki/Fisher's_method
	 * 
	 * @return
	 */
	public double pValueFisherChi2() {
		if (size() <= 0) return 1.0;

		double sum = 0.0;
		int count = 0;
		for (int i = 0; i < size(); i++) {
			double pvalue = getPvalue(i);

			// If pvalue == 0, it produces an error (log will be -Inf)
			if (pvalue > 0) {
				sum += Math.log(pvalue);
				count++;
			}
		}

		// Nothing added? 
		if (count <= 0) return 1.0;

		// Get new p-value
		double chi2 = -2.0 * sum;
		int k = 2 * count;
		double pValue = chiSquareCDFComplementary(chi2, k); // 1 - ChiSquareCDF_{k}( chi2 )

		return pValue;
	}

	/**
	 * Get minimum pvalue
	 * @return
	 */
	public double pValueMin() {
		double min = 1.0;
		for (int i = 0; i < size(); i++)
			min = Math.min(min, getPvalue(i));
		return min;
	}

	/**
	 * Combine p-values using Simes's procedure
	 * 
	 * References: http://biomet.oxfordjournals.org/content/73/3/751
	 * 
	 * @return A combined p-value
	 */
	public double pValueSimes() {
		if (size() <= 0) return 1.0;

		// Count non-zero p-values (we treat zero p-values as errors, so we skip them) 
		int total = 0;
		for (int i = 0; i < size(); i++)
			if (getPvalue(i) > 0) total++;
		double tot = total;

		// No p-value? => We are done
		if (total <= 0) return 1.0;

		sort(); // Sort collection

		// Perform Simes's method
		int count = 0;
		double pSimesMin = 1.0;
		for (int i = 0; i < size(); i++) {
			double pvalue = getPvalue(i);

			// We treat zero p-values as errors, so we skip them 
			if (pvalue > 0) {
				count++;
				double pSimes = (pvalue * tot) / count;
				pSimesMin = Math.min(pSimesMin, pSimes);
			}
		}

		return pSimesMin;
	}

	/**
	 * Combine p-values using Stouffer's Z-score method
	 * 
	 * References: http://en.wikipedia.org/wiki/Fisher's_method  (scroll down to Stouffer's method)
	 * 
	 * @return A combined p-value
	 */
	public double pValueZScore() {
		if (size() <= 0) return 1.0;

		double sum = 0.0;
		int count = 0;
		for (int i = 0; i < size(); i++) {
			double pvalue = getPvalue(i);

			// If pvalue == 0, it produces an error (normal inverse is -Inf)
			if (pvalue > 0) {
				double z = Stat.normalInverseCDF(pvalue);
				sum += z;

				count++;
			}
		}

		// Nothing added? 
		if (count <= 0) return 1.0;

		// Get new p-value
		double zsum = sum / Math.sqrt(count);
		double pValue = Stat.normalCDF(0.0, 1.0, zsum);

		return pValue;
	}

	/**
	 * Get pvalue quantile
	 * @return
	 */
	public double quantile(double quantile) {
		if (quantile < 0.0 || quantile > 1.0) throw new RuntimeException("Quantile out of range: " + quantile + " .Expected range [0.0 , 1.0].");
		if (size() <= 0) return 1.0;

		sort(); // Sort collection
		int num = (int) (quantile * size());
		return pValues.get(num);
	}

	public void setGeneId(String geneId) {
		this.geneId = geneId;
	}

	public int size() {
		return pValues.size();
	}

	void sort() {
		if (!sorted) {
			pValues.sort();
			sorted = true;
		}
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();

		sb.append(geneId + ":\t");
		for (int i = 0; i < size(); i++)
			sb.append(String.format(" %.2e", getPvalue(i)));

		return sb.toString();
	}

}
