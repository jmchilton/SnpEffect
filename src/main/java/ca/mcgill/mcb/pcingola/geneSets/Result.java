package ca.mcgill.mcb.pcingola.geneSets;

import java.util.LinkedList;
import java.util.List;

import org.apfloat.Apcomplex;
import org.apfloat.Apfloat;

/**
 * Store a result form a greedy search algorithm
 *
 * @author Pablo Cingolani
 */
public class Result implements Comparable<Result> {

	private List<GeneSet> list;
	private Apfloat pValue;
	private double geneSetCount;

	public Result() {
		list = null;
		pValue = Apcomplex.ONE;
	}

	public Result(GeneSet geneSet, Apfloat pValue, double geneSetCount) {
		list = new LinkedList<GeneSet>();
		list.add(geneSet);
		this.pValue = pValue;
		this.geneSetCount = geneSetCount;
	}

	public Result(List<GeneSet> list, double pValue, double geneSetCount) {
		this.list = list;
		this.pValue = new Apfloat(pValue);
		this.geneSetCount = geneSetCount;
	}

	@Override
	public int compareTo(Result res) {
		if ((pValue != null) && (res.pValue != null)) return pValue.compareTo(res.pValue);
		return 0;
	}

	public double getGeneSetCount() {
		return geneSetCount;
	}

	public GeneSet getLatestGeneSet() {
		if ((list == null) || (list.size() == 0)) { return null; }
		return list.get(list.size() - 1);
	}

	public List<GeneSet> getList() {
		return list;
	}

	public Apfloat getPvalue() {
		return pValue;
	}

	/**
	 * P-Value adjusted using FDR
	 * @return
	 */
	public double getPvalueAdjusted() {
		double adj = 1.0;
		for (int i = 0; i < list.size(); i++)
			adj *= geneSetCount / (i + 1);
		return Math.min(1.0, pValue.doubleValue() * adj);
	}

	public double getPvalueDouble() {
		return pValue.doubleValue();
	}

	public void setGeneSetCount(double geneSetCount) {
		this.geneSetCount = geneSetCount;
	}

	public void setList(List<GeneSet> list) {
		this.list = list;
	}

	public void setPvalue(Apfloat pValue) {
		this.pValue = pValue;
	}

	public void setPvalueDouble(double pValue) {
		this.pValue = new Apfloat(pValue);
	}
}
