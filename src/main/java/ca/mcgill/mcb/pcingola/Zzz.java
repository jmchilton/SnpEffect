package ca.mcgill.mcb.pcingola;

import java.util.Random;

import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Genome;
import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.SnpEffectPredictor;
import ca.mcgill.mcb.pcingola.util.GprSeq;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * Simple test program
 * @author pcingola
 */
public class Zzz {

	public static void main(String[] args) {
		Config config = new Config("testHg3765Chr22");
		Timer.show("Loading predictor");
		SnpEffectPredictor snpEffectPredictor = config.loadSnpEffectPredictor();
		Timer.show("Done");

		Random random = new Random(20130214);

		Genome genome = snpEffectPredictor.getGenome();
		for (Gene g : genome.getGenes()) {
			if (g.isProteinCoding()) {
				System.out.println(g.getGeneName());

				for (Transcript t : g) {
					System.out.println("\t" + t.getId());

					for (int i = t.getStart(); i < t.getEnd(); i++) {
						// Create a fake SNP. Random REF and ALT bases
						char alt, ref = GprSeq.randBase(random);
						do {
							alt = GprSeq.randBase(random);
						} while (ref == alt);

						System.out.println("\t\t" + i + "\t" + ref + "\t" + alt);

						String refStr = ref + "";
						String altStr = alt + "";
						SeqChange sc = new SeqChange(t.getChromosome(), i, refStr, altStr, 1, "", -1, -1);

						Transcript tnew = t.apply(sc);
					}
				}
			}
		}
	}
}
