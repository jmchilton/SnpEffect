package ca.mcgill.mcb.pcingola.gsa;

import gnu.trove.list.array.TDoubleArrayList;

/**
 * A list of pvalues for a gene
 * 
 * @author pcingola
 */
public class GenePvalueList {

	public enum PvalueSummary {
		MIN, AVG, AVG10
	};

	String geneId;
	TDoubleArrayList pValues;

	public GenePvalueList() {
		pValues = new TDoubleArrayList();
	}

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
