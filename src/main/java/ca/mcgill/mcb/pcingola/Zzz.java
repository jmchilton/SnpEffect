package ca.mcgill.mcb.pcingola;

import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Genome;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * Simple test program
 * @author pcingola
 */
public class Zzz {

	public static void main(String[] args) {
		if (args.length != 1) {
			System.err.println("Usage: " + Zzz.class.getSimpleName() + " genomeVer");
			System.exit(-1);
		}
		String genomeVer = args[0];

		Timer.showStdErr("Loading database " + genomeVer);

		Config config = new Config(genomeVer);
		config.loadSnpEffectPredictor();
		Genome genome = config.getGenome();

		Timer.showStdErr("Calculate");
		System.out.println("gene.id\ttr.id\te.rank\tg.numChilds\ttr.numChilds\tg.size\ttr.size\te.size");
		for (Gene g : genome.getGenes()) {
			if (g.isProteinCoding()) {
				for (Transcript tr : g) {
					if (tr.isProteinCoding()) {
						for (Exon e : tr) {
							System.out.println(g.getId() + "\t" + tr.getId() + "\t" + e.getRank() + "\t" + g.numChilds() + "\t" + tr.numChilds() + "\t" + g.size() + "\t" + tr.size() + "\t" + e.size());
						}
					}
				}
			}
		}

		Timer.showStdErr("Done");
	}
}
