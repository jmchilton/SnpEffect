package ca.mcgill.mcb.pcingola;

import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Genome;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.stats.IntStats;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * Simple test program
 * @author pcingola
 */
public class Zzz {

	public static void main(String[] args) {
		String genomeVer = "testHg3766Chr1";

		Timer.showStdErr("Loading database");

		Config config = new Config(genomeVer);
		config.loadSnpEffectPredictor();
		Genome genome = config.getGenome();

		Timer.showStdErr("Calculate");
		IntStats numTr = new IntStats();
		IntStats numExon = new IntStats();
		IntStats trSize = new IntStats();
		IntStats exonSize = new IntStats();

		for (Gene g : genome.getGenes()) {
			Gpr.debug(g.getId());
			numTr.sample(g.numChilds());

			for (Transcript tr : g) {
				trSize.sample(tr.size());
				numExon.sample(tr.numChilds());
				for (Exon e : tr)
					exonSize.sample(e.size());
			}
		}

		System.err.println("Number of genes: \n" + genome.getGenes().size());
		System.err.println("Number of transcripts: \n" + numTr);
		System.err.println("Number of exons: \n" + numExon);
		System.err.println("Transcript size: \n" + trSize);
		System.err.println("Exon sizes: \n" + exonSize);
		Timer.showStdErr("Done");
	}
}
