package ca.mcgill.mcb.pcingola;

import java.util.Random;

import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.interval.Upstream;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.SnpEffectPredictor;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;

public class Zzz {

	public static final int READ_LENGTH = 50;
	public static final int NUM_READS_MARKER = 5;

	Config config;
	SnpEffectPredictor sep;
	Random random;
	StringBuilder out;

	String outFile = Gpr.HOME + "/fly_pvuseq/rand_up.bed";

	public static void main(String[] args) {

		String fileName = "adsf/adsfadf*\' 4576467  > 0";
		System.out.println(fileName);
		System.out.println(Gpr.sanityzeFileName(fileName));

		//		Zzz zzz = new Zzz();
		//		zzz.run();
	}

	public Zzz() {
		Timer.showStdErr("Loading");
		config = new Config("BDGP5.69");
		sep = config.loadSnpEffectPredictor();
		Timer.showStdErr("Building");
		sep.buildForest();
		Timer.showStdErr("Done");
	}

	/**
	 * Random read on this marker
	 * @param m
	 */
	void randReads(Marker m) {
		randReads(m.getChromosomeName(), m.getStart(), m.getEnd(), m.getStrand(), m.idChain());
	}

	void randReads(String chrName, int mstart, int mend, int strand, String id) {
		int size = mend - mstart;
		if (size <= 0) return;

		for (int i = 0; i < NUM_READS_MARKER; i++) {
			// Up
			int r1 = random.nextInt(size);
			int r2 = random.nextInt(size);

			int rand;
			if (strand >= 0) rand = Math.max(r1, r2);
			else rand = Math.min(r1, r2);

			int start = mstart + rand;
			int end = start + READ_LENGTH;

			out.append(chrName + "\t" + start + "\t" + end + "\t" + id + "\n");
		}
	}

	/**
	 * Random read on every exon
	 * @param tr
	 */
	void randReadsExons(Transcript tr) {
		tr.rankExons();

		for (Exon e : tr.sortedStrand()) {
			if (tr.getStrand() != e.getStrand()) {
				System.err.println("WTF!?!?\n\t" + tr.getId() + "\t" + tr.getStrand() + "\n\t" + e.getId() + "\t" + e.getStrand());
				continue;
			}

			for (int i = 0; i < NUM_READS_MARKER; i++)
				randReads(e);
		}
	}

	/**
	 * Read upstream & downstream
	 * @param tr
	 */
	void randReadsUpDownStream(Gene g) {
		// Pick ppstream region that is upstream of every transcript
		int start = Integer.MAX_VALUE, end = Integer.MAX_VALUE;
		Transcript t = null;
		for (Transcript tr : g) {
			Upstream up = tr.getUpstream();
			if (up != null) {
				start = Math.min(start, up.getStart());
				end = Math.min(end, up.getEnd());
				t = tr;
			}
		}

		// if (start < end) randReads(g.getChromosomeName(), start, end, g.getStrand(), t.getUpstream().idChain());
		if (start < end) unifReads(g.getChromosomeName(), start, end, g.getStrand(), t.getUpstream().idChain());
		if (start == end) Gpr.debug(start + "\t" + end);
	}

	// Save
	void run() {
		random = new Random(20130601);
		out = new StringBuilder();

		for (Gene g : sep.getGenome().getGenes()) {
			//if (g.isStrandPlus()) randReadsUpDownStream(g);
			for (Transcript tr : g)
				randReadsExons(tr);
		}

		Timer.showStdErr("Saving to " + outFile);
		Gpr.toFile(outFile, out);
	}

	/**
	 * Uniform 'read' covering the region
	 * @param chrName
	 * @param mstart
	 * @param mend
	 * @param strand
	 * @param id
	 */
	void unifReads(String chrName, int mstart, int mend, int strand, String id) {
		int size = mend - mstart;

		int step = size / NUM_READS_MARKER;
		if (step <= 0) step = 1;

		for (int pos = mstart; pos <= mend; pos += step) {
			int start = pos;
			int end = start + READ_LENGTH;
			out.append(chrName + "\t" + start + "\t" + end + "\t" + id + "\n");
		}
	}
}
