package ca.mcgill.mcb.pcingola.gsa;

import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;

import java.util.ArrayList;

import ca.mcgill.mcb.pcingola.interval.Chromosome;

/**
 * A list of <chromosome, position, pvalues> 
 * 
 * @author pcingola
 */
public class ChrPosPvalueList {

	ArrayList<Chromosome> chromosomes;
	TIntArrayList starts;
	TIntArrayList ends;
	TDoubleArrayList pValues;

	public ChrPosPvalueList() {
		chromosomes = new ArrayList<Chromosome>();
		starts = new TIntArrayList();
		ends = new TIntArrayList();
		pValues = new TDoubleArrayList();
	}

	public void add(Chromosome chr, int start, int end, double pvalue) {
		chromosomes.add(chr);
		starts.add(start);
		ends.add(start);
		pValues.add(pvalue);
	}

	public Chromosome getChromosome(int index) {
		return chromosomes.get(index);
	}

	public String getChromosomeName(int index) {
		return chromosomes.get(index).getId();
	}

	public int getEnd(int index) {
		return ends.get(index);
	}

	public double getPvalue(int index) {
		return pValues.get(index);
	}

	public int getStart(int index) {
		return starts.get(index);
	}

	public int size() {
		return chromosomes.size();
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();

		for (int i = 0; i < size(); i++)
			sb.append(getChromosomeName(i) + "\t" + getStart(i) + "\t" + getPvalue(i) + "\n");

		return sb.toString();
	}

}
