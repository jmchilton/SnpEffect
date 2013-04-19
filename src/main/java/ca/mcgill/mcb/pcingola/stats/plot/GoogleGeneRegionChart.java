package ca.mcgill.mcb.pcingola.stats.plot;

import java.util.ArrayList;

import ca.mcgill.mcb.pcingola.interval.Downstream;
import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Intron;
import ca.mcgill.mcb.pcingola.interval.Upstream;
import ca.mcgill.mcb.pcingola.interval.Utr3prime;
import ca.mcgill.mcb.pcingola.interval.Utr5prime;
import ca.mcgill.mcb.pcingola.stats.CoverageByType;
import ca.mcgill.mcb.pcingola.stats.PosStats;
import ca.mcgill.mcb.pcingola.util.Gpr;

public class GoogleGeneRegionChart {

	CoverageByType coverageByType;
	GoogleLineChart lineChart;
	String name;
	String[] types = { Upstream.class.getSimpleName(), Utr5prime.class.getSimpleName(), Exon.class.getSimpleName(), Intron.class.getSimpleName(), Utr3prime.class.getSimpleName(), Downstream.class.getSimpleName() };

	public GoogleGeneRegionChart(CoverageByType coverageByType, String name) {
		this.coverageByType = coverageByType;
		this.name = name;
		init();
	}

	int addCol(int nullBefore, String type) {
		ArrayList<String> values = createCol(nullBefore, type);
		lineChart.addColumn(type, values);
		return values.size();
	}

	ArrayList<String> createCol(int nullBefore, String type) {
		ArrayList<String> col = new ArrayList<String>();
		PosStats posStats = coverageByType.get(type);
		Gpr.debug(type + " Size: " + posStats.size());

		for (int i = 0; i < nullBefore; i++)
			col.add(null);

		for (int i = 0; i < posStats.size(); i++)
			col.add("" + posStats.getCount(i));

		Gpr.debug(col.size());
		return col;
	}

	void init() {
		lineChart = new GoogleLineChart(name);
	}

	public String toStringHtmlBody() {
		int nullBefore = 0;
		for (String type : types)
			nullBefore = addCol(nullBefore, type);

		return lineChart.toStringHtmlHeader() + lineChart.toStringHtmlBody();
	}

	public String toStringHtmlHeader() {
		StringBuilder sb = new StringBuilder();
		return sb.toString();
	}

}
