package ca.mcgill.mcb.pcingola;

import java.util.Random;

import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.SnpEffectPredictor;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;

public class Zzz {

	public static final int READ_LENGTH = 50;
	public static final int NUM_READS_EXON = 5;

	Config config;
	SnpEffectPredictor sep;

	public static void main(String[] args) {
		Zzz zzz = new Zzz();
		zzz.run();
	}

	public Zzz() {
		Timer.showStdErr("Loading");
		config = new Config("BDGP5.69");
		sep = config.loadSnpEffectPredictor();
		Timer.showStdErr("Done");
	}

	void run() {
		Random random = new Random(20130601);
		StringBuilder out = new StringBuilder();

		for (Gene g : sep.getGenome().getGenes()) {
			for (Transcript tr : g) {
				tr.rankExons();

				for (Exon e : tr.sortedStrand()) {

					if (tr.getStrand() != e.getStrand()) {
						System.err.println("WTF!?!?\n\t" + tr.getId() + "\t" + tr.getStrand() + "\n\t" + e.getId() + "\t" + e.getStrand());
						continue;
					}

					for (int i = 0; i < NUM_READS_EXON; i++) {
						// Up
						int r1 = random.nextInt(e.size());
						int r2 = random.nextInt(e.size());

						int rand;
						if (e.isStrandPlus()) rand = Math.max(r1, r2);
						else rand = Math.min(r1, r2);

						int start = e.getStart() + rand;
						int end = start + READ_LENGTH;
						out.append(e.getChromosomeName() + "\t" + start + "\t" + end + "\n");
					}
				}
			}
		}

		// Save
		String outFile = Gpr.HOME + "/fly_pvuseq/rand_up.bed";
		Timer.showStdErr("Saving to " + outFile);
		Gpr.toFile(outFile, out);
	}
}
