package ca.mcgill.mcb.pcingola;

import ca.mcgill.mcb.pcingola.codons.CodonTable;
import ca.mcgill.mcb.pcingola.codons.CodonTables;
import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.SnpEffectPredictor;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;

public class Zzz {

	Config config;
	SnpEffectPredictor sep;
	byte bases[];

	public static void main(String[] args) {

		String cds = "GATATCTGCTTGGTATCTTCTCAATATCTTGACCATCTGTGACAATTTTAATCCTCATTT";
		CodonTable ct = CodonTables.getInstance().getTable(CodonTables.STANDARD_TABLE_NAME);

		for (int i = 0; i < cds.length(); i++) {
			String c = cds.substring(i);
			String aa = ct.aa(c);
			System.out.println("AA: " + aa + "\t" + c);
		}

		//		Zzz zzz = new Zzz();
		//		zzz.load("testHg3770Chr22");
		//		zzz.run();
	}

	public Zzz() {
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

	void set(Gene g, int pos) {
		int i = pos - g.getStart();
		bases[i] = 1;
	}
}
