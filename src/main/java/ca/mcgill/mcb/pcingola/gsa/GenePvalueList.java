package ca.mcgill.mcb.pcingola.gsa;

import flanagan.analysis.Stat;
import gnu.trove.list.array.TDoubleArrayList;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * A list of pvalues for a gene
 * 
 * @author pcingola
 */
public class GenePvalueList {

	public enum PvalueSummary {
		MIN, AVG, AVG10, FISHER_CHI_SQUARE
	};

	String geneId;
	TDoubleArrayList pValues;

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

	public GenePvalueList() {
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

		pValues.sort(); // Sort collection

		// Sum smallest values
		double sum = 0.0;
		int max = Math.min(size(), topN);
		for (int i = 0; i < max; i++)
			sum += getPvalue(i);

		return sum / max;
	}

	/**
	 * Combine p-values using Fisher's method
	 * 
	 * References: http://en.wikipedia.org/wiki/Fisher's_method
	 * 
	 * @return
	 */
	double pValueFisherChi2() {
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
		Gpr.debug("Chi2 : " + chi2 + "\tk: " + k + "\tpval: " + pValue);

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

	public void setGeneId(String geneId) {
		this.geneId = geneId;
	}

	public int size() {
		return pValues.size();
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
