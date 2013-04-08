package ca.mcgill.mcb.pcingola;

import ca.mcgill.mcb.pcingola.interval.Cds;
import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.SnpEffectPredictor;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * Tast class
 * @author pablocingolani
 */
public class Zzz {
	int exonsNonZeroFrame = 0;

	/**
	 * Main
	 * @param args
	 */
	public static void main(String[] args) {
		// String gen = "testHg3763Chr20";
		String gen = "hg19";
		Zzz zzz = new Zzz();
		zzz.run(gen);
	}

	void run(String gen) {
		// Load data
		Timer.showStdErr("Loading..");
		Config config = new Config(gen);
		SnpEffectPredictor snpEffectPredictor = config.loadSnpEffectPredictor();

		// Analyze genes
		int countErrsTr = 0, totalTr = 0;
		for (Gene g : snpEffectPredictor.getGenome().getGenes()) {
			for (Transcript tr : g) {
				if (!checkExonFrames(tr)) countErrsTr++;
				totalTr++;
			}
		}

		Timer.showStdErr("\n\tErrors                    : " + countErrsTr // 
				+ "\n\tTranscripts               : " + totalTr //
				+ "\n\tExons with non-zero frame : " + exonsNonZeroFrame //
		);
	}

	/**
	 * Check that exon frames match
	 * @param tr
	 * @return
	 */
	boolean checkExonFrames(Transcript tr) {
		boolean ok = true;

		int bases = 0;
		if (tr.isStrandPlus()) {
			Cds codingPart = new Cds(tr, tr.getCdsStart(), tr.getCdsEnd(), tr.getStrand(), tr.getId());
			for (Exon e : tr.sortedStrand()) {
				if (e.intersects(codingPart)) {

					// Calculate frame
					int frame = (bases % 3);
					//if (frame > 0) frame = 3 - frame;  // USE THIS LINE FOR ENSEMBL, BUT NOT FOR UCSC REFSEQ

					if (e.getFrame() > 0) {
						exonsNonZeroFrame++;

						// Does it match?
						if (frame != e.getFrame()) {
							System.out.println("Frame error: " + tr.getId() + ", exon " + e.getRank() + "\tFrame (ann): " + e.getFrame() + "\tFrame (calc) : " + frame);
							ok = false;
						}
					}

					// Update base count
					bases += e.getEnd() - Math.max(codingPart.getStart(), e.getStart()) + 1;
				}
			}
		}

		return ok;
	}
}
