package ca.mcgill.mcb.pcingola;

import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.SnpEffectPredictor;
import ca.mcgill.mcb.pcingola.util.Timer;

public class Zzz {

	public static final int READ_LENGTH = 50;
	public static final int NUM_READS_MARKER = 5;

	Config config;
	SnpEffectPredictor sep;

	public static void main(String[] args) {

		String str = "1\t2\t\t\t\t3";
		String s[] = str.split("\t");
		for (int i = 0; i < s.length; i++)
			System.out.println(i + "\t" + s[i]);

		//		String genome = "hg19";
		//		Zzz zzz = new Zzz(genome);
		//		zzz.run();
	}

	public Zzz(String genome) {
		config = new Config(genome);
	}

	public void run() {
		Timer.showStdErr("Loading");
		sep = config.loadSnpEffectPredictor();
		Timer.showStdErr("Building");
		sep.buildForest();
		Timer.showStdErr("Done");

		for (Gene gene : sep.getGenome().getGenes()) {
			for (Transcript tr : gene) {
				if ((tr.numChilds() > 100) || (gene.numChilds() > 100)) System.out.println(tr.numChilds() + "\t" + gene.numChilds() + "\t" + tr.getId() + "\t" + gene.getId() + "\t" + gene.getGeneName());
			}
		}
	}
}
