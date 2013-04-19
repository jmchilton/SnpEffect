package ca.mcgill.mcb.pcingola.snpEffect.commandLine;

import java.util.ArrayList;

import ca.mcgill.mcb.pcingola.interval.Downstream;
import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Intron;
import ca.mcgill.mcb.pcingola.interval.Upstream;
import ca.mcgill.mcb.pcingola.interval.Utr3prime;
import ca.mcgill.mcb.pcingola.interval.Utr5prime;
import ca.mcgill.mcb.pcingola.stats.CoverageByType;
import ca.mcgill.mcb.pcingola.stats.plot.GoogleGeneRegionChart;
import ca.mcgill.mcb.pcingola.stats.plot.GoogleLineChart;
import ca.mcgill.mcb.pcingola.util.Gpr;

public class Zzz {

	public static void main(String[] args) {
		Zzz zzz = new Zzz();
		zzz.run();
	}

	void run() {
		CoverageByType coverageByType = new CoverageByType();
		coverageByType.getOrCreate(Upstream.class.getSimpleName()).rand(100, 100);
		coverageByType.getOrCreate(Utr5prime.class.getSimpleName()).rand(100, 100);
		coverageByType.getOrCreate(Exon.class.getSimpleName()).rand(100, 100);
		coverageByType.getOrCreate(Intron.class.getSimpleName()).rand(100, 100);
		coverageByType.getOrCreate(Utr3prime.class.getSimpleName()).rand(100, 100);
		coverageByType.getOrCreate(Downstream.class.getSimpleName()).rand(100, 100);

		GoogleGeneRegionChart grc = new GoogleGeneRegionChart(coverageByType, "Gene");
		Gpr.toFile(Gpr.HOME + "/z.html", grc.toStringHtmlHeader() + grc.toStringHtmlBody());
		Gpr.debug("Done!");
	}

	void runLineChart() {
		GoogleLineChart glc = new GoogleLineChart("Name");
		glc.setTitle("TITLE!");
		glc.setWidth(1200);
		glc.setvAxis("Y-axis label");
		glc.sethAxis("X-axis label");

		for (int j = 0; j < 5; j++) {
			ArrayList<String> col = new ArrayList<String>();
			for (int i = 0; i < 100; i++)
				col.add(Math.random() + "");
			glc.addColumn("Col_Number" + j, col);
		}

		Gpr.toFile(Gpr.HOME + "/z.html", glc.toStringHtmlHeader() + glc.toStringHtmlBody());
		Gpr.debug("Done!");
	}
}
