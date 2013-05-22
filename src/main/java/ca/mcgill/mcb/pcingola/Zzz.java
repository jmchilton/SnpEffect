package ca.mcgill.mcb.pcingola;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;

import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.SnpEffectPredictor;
import ca.mcgill.mcb.pcingola.stats.CoverageByType;
import ca.mcgill.mcb.pcingola.stats.PosStats;
import ca.mcgill.mcb.pcingola.stats.plot.GoogleLineChart;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;

public class Zzz {

	Config config;
	SnpEffectPredictor sep;
	byte bases[];

	public static void main(String[] args) {
		Zzz zzz = new Zzz();
		zzz.runTrBugs();
	}

	public Zzz() {
	}

	GoogleLineChart lineChart(String chr, PosStats posstats) {
		GoogleLineChart glc = new GoogleLineChart("Histogram");
		glc.setTitle("Chromosome: " + chr);
		glc.setWidth(1200);
		glc.setvAxis("Transcripts");
		glc.sethAxis("Position");

		ArrayList<String> col = new ArrayList<String>();
		for (int i = 0; i < posstats.size(); i++)
			col.add(posstats.getCount(i) + "");
		glc.addColumn("Chr" + chr, col);

		return glc;
	}

	void load(String genVer) {
		Timer.showStdErr("Loading");
		config = new Config(genVer);
		sep = config.loadSnpEffectPredictor();
		System.out.println(sep.getGenome());
	}

	void run() {
		Timer.showStdErr("Checking");

		for (Gene g : sep.getGenome().getGenes()) {
			// Initialize
			bases = new byte[g.size()];
			for (int i = 0; i < bases.length; i++)
				bases[i] = 0;

			run(g);

			// Count
			int count = 0;
			for (int i = 0; i < bases.length; i++)
				if (bases[i] > 0) count++;

			System.out.println(g.getGeneName() + "\t" + count);
		}

		Timer.showStdErr("Done");
	}

	void run(Gene g) {
		for (Transcript tr : g) {
			int cdsStart, cdsEnd;

			if (tr.isStrandMinus()) {
				cdsStart = tr.getCdsEnd();
				cdsEnd = tr.getCdsStart();
			} else {
				cdsStart = tr.getCdsStart();
				cdsEnd = tr.getCdsEnd();
			}

			Gpr.debug(tr.getStrand() + "\t" + cdsStart + "\t" + cdsEnd);

			if (tr.isProteinCoding()) {
				for (Exon e : tr) {
					// Set bases in coding part
					int min = Math.max(cdsStart, e.getStart());
					int max = Math.min(cdsEnd, e.getEnd());

					for (int i = min; i <= max; i++)
						set(g, i);

					// Set splice sites
				}
			}
		}
	}

	public void runTrBugs() {
		//		String genome = "testHg3766Chr1";
		String genome = "GRCh37.70";
		Timer.showStdErr("Loading");

		Config config = new Config(genome);
		SnpEffectPredictor sep = config.loadSnpEffectPredictor();

		CoverageByType coverageByType = new CoverageByType();

		StringBuilder sb = new StringBuilder();
		int count = 1;
		for (Gene g : sep.getGenome().getGenes()) {
			for (Transcript tr : g) {
				if (tr.isProteinCoding()) {
					boolean hasError = false;

					if (tr.isErrorProteinLength()) {
						hasError = true;
					}

					if (tr.isErrorStopCodonsInCds()) {
						hasError = true;
					}

					if (tr.isErrorStopCodon()) {
						// hasError = true;
					}

					if (tr.isErrorStartCodon()) {
						hasError = true;
					}

					if (hasError) {
						int perc = (int) (100.0 * (((double) tr.getStart()) / tr.getChromosome().getEnd()));

						String out = count //
								+ "\t" + g.getGeneName() //
								+ "\t" + g.getId() //
								+ "\t" + tr.getId() //
								+ "\t" + tr.getChromosomeName() //
								+ "\t" + tr.getStart() //
								+ "\t" + tr.getEnd() //
								+ "\t" + perc//
						;

						coverageByType.getOrCreate(tr.getChromosomeName()).sample(tr, tr.getChromosome());

						System.out.println(out);
						sb.append(out + "\n");

						count++;
					}
				}
			}
		}

		//---
		// Save data to file
		//---
		String fileName = Gpr.HOME + "/tr_bugs.txt";
		Timer.showStdErr("Saving to file " + fileName);
		Gpr.toFile(fileName, sb.toString());

		//---
		// Show charts
		//---
		ArrayList<String> chrs = new ArrayList<String>();
		chrs.addAll(coverageByType.keySet());
		Collections.sort(chrs);
		LinkedList<GoogleLineChart> glcs = new LinkedList<GoogleLineChart>();
		for (String chr : chrs) {
			if (chr.length() <= 2) {
				PosStats posstats = coverageByType.get(chr);
				glcs.add(lineChart(chr, posstats));
			}
		}

		// Save charts to HTML file
		StringBuilder html = new StringBuilder();
		for (GoogleLineChart glc : glcs)
			html.append(glc.toStringHtmlHeader());
		for (GoogleLineChart glc : glcs)
			html.append(glc.toStringHtmlBody());
		fileName = Gpr.HOME + "/tr_bugs.html";
		Timer.showStdErr("Saving to file " + fileName);
		Gpr.toFile(fileName, html.toString());

		Timer.showStdErr("Done.");
	}

	void set(Gene g, int pos) {
		int i = pos - g.getStart();
		bases[i] = 1;
	}
}
