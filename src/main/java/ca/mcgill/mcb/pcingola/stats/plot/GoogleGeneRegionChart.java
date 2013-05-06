package ca.mcgill.mcb.pcingola.stats.plot;

import java.util.ArrayList;

import ca.mcgill.mcb.pcingola.interval.Downstream;
import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Intergenic;
import ca.mcgill.mcb.pcingola.interval.Intron;
import ca.mcgill.mcb.pcingola.interval.Upstream;
import ca.mcgill.mcb.pcingola.interval.Utr3prime;
import ca.mcgill.mcb.pcingola.interval.Utr5prime;
import ca.mcgill.mcb.pcingola.stats.CoverageByType;
import ca.mcgill.mcb.pcingola.stats.PosStats;

public class GoogleGeneRegionChart {

	String[] types = { Intergenic.class.getSimpleName() //
			, Upstream.class.getSimpleName() //
			, Utr5prime.class.getSimpleName() //
			, Exon.class.getSimpleName() //
			, Intron.class.getSimpleName() //
			, Utr3prime.class.getSimpleName() //
			, Downstream.class.getSimpleName() //
	};

	CoverageByType coverageByType;
	GoogleLineChart lineChart;
	String name, header, body;

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

	void createChart() {
		int nullBefore = 0;
		for (String type : types)
			nullBefore = addCol(nullBefore, type);

		header = lineChart.toStringHtmlHeader();
		body = lineChart.toStringHtmlBody();
	}

	ArrayList<String> createCol(int nullBefore, String type) {
		ArrayList<String> col = new ArrayList<String>();
		PosStats posStats = coverageByType.get(type);

		for (int i = 0; i < nullBefore; i++)
			col.add(null);

		for (int i = 0; i < posStats.size(); i++)
			col.add("" + posStats.getCount(i));

		return col;
	}

	void init() {
		lineChart = new GoogleLineChart(name);
	}

	public String toStringHtmlBody() {
		if (body == null) createChart();
		return body;
	}

	public String toStringHtmlHeader() {
		if (header == null) createChart();
		return header;
	}

}
