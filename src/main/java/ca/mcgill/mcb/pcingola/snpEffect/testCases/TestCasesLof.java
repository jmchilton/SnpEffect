package ca.mcgill.mcb.pcingola.snpEffect.testCases;

import java.util.LinkedList;

import junit.framework.Assert;
import junit.framework.TestCase;
import net.sf.samtools.util.RuntimeEOFException;
import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Markers;
import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.interval.SpliceSite;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.LossOfFunction;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * Test Loss of Function prediction
 * 
 * @author pcingola
 */
public class TestCasesLof extends TestCase {

	public static boolean debug = true;

	Config config;

	public TestCasesLof() {
		super();
	}

	/**
	 * Create change effects
	 * @param seqChange
	 * @param effectType
	 * @param exon
	 * @return
	 */
	LinkedList<ChangeEffect> changeEffects(SeqChange seqChange, EffectType effectType, Marker marker) {
		ChangeEffect changeEffect = new ChangeEffect(seqChange);
		changeEffect.set(marker, effectType, "");
		LinkedList<ChangeEffect> changeEffects = new LinkedList<ChangeEffect>();
		changeEffects.add(changeEffect);
		return changeEffects;

	}

	/**
	 * Check that LOF works for a given transcript
	 * @param gene
	 * @param tr
	 */
	void checkLof(Transcript tr) {
		// checkLofSplice(tr);
		//		checkLofStartLost(tr);
		checkLofExonDeleted(tr);

		//		checkLofFrameShift(tr);
	}

	/**
	 * Check that exon deleted is LOF
	 * @param tr
	 */
	void checkLofExonDeleted(Transcript tr) {
		for (Exon ex : tr) {
			SeqChange seqChange = new SeqChange(tr.getChromosome(), ex.getStart(), ex.getEnd(), ""); // Create a seqChange
			Gpr.debug("SeqChange: " + seqChange);
			LinkedList<ChangeEffect> changeEffects = changeEffects(seqChange, EffectType.EXON_DELETED, ex);
			LossOfFunction lof = new LossOfFunction(changeEffects);
			boolean islof = lof.isLof();
			Assert.assertEquals(true, islof);
		}
	}

	void checkLofFrameShift(Transcript tr) {
		throw new RuntimeException("Unimplemented!");
	}

	void checkLofSplice(Transcript tr) {
		// All transcripts in exon
		int exonNum = 0;
		for (Exon ex : tr.sortedStrand()) {
			checkSpliceDonor(tr, ex, exonNum);
			checkSpliceAcceptor(tr, ex, exonNum);
			exonNum++;
		}
	}

	/**
	 * Check that START_LOST is LOF
	 * @param tr
	 */
	void checkLofStartLost(Transcript tr) {
		// Find start codon position
		int pos = tr.getCdsStart();
		SeqChange seqChange = new SeqChange(tr.getChromosome(), pos, pos, ""); // Create a seqChange

		// Finr exon
		Exon exon = null;
		for (Exon ex : tr)
			if (ex.intersects(pos)) exon = ex;
		if (exon == null) throw new RuntimeEOFException("Cannot find first exon for transcript " + tr.getId());

		// Create a LOF object and analyze the effect
		LinkedList<ChangeEffect> changeEffects = changeEffects(seqChange, EffectType.START_LOST, exon);
		LossOfFunction lof = new LossOfFunction(changeEffects);
		boolean islof = lof.isLof();
		Assert.assertEquals(true, islof);
	}

	/**
	 * Check that Core Splice Site acceptors are considered LOF
	 * @param tr
	 * @param ex
	 * @param exonNum
	 */
	void checkSpliceAcceptor(Transcript tr, Exon ex, int exonNum) {
		int step = tr.isStrandPlus() ? -1 : +1;
		int intronNum = exonNum - 1; // We care about the intron before

		if (ex.getRank() > 1) {
			// Position
			int posDonor = tr.isStrandPlus() ? ex.getStart() : ex.getEnd();
			posDonor += step;

			// Splice site size
			int maxSize = Math.min(tr.intronSize(intronNum), SpliceSite.CORE_SPLICE_SITE_SIZE);
			if (debug) Gpr.debug("Intron size: " + tr.intronSize(intronNum));
			if (maxSize <= 0) throw new RuntimeEOFException("Max splice size is " + maxSize);

			//---
			// For all position on splice site donor positions, make sure it is LOF
			//---
			for (int pos = posDonor, i = 0; i < maxSize; i++, pos += step) {
				SeqChange seqChange = new SeqChange(tr.getChromosome(), pos, pos, ""); // Create a seqChange
				Marker marker = findMarker(seqChange, EffectType.SPLICE_SITE_ACCEPTOR, null, ex);
				LinkedList<ChangeEffect> changeEffects = changeEffects(seqChange, EffectType.SPLICE_SITE_ACCEPTOR, marker); // Create a SPLICE_SITE_ACCEPTOR effect
				if (debug) Gpr.debug("SeqChange:" + seqChange);

				// Create a LOF object and analyze the effect
				LossOfFunction lof = new LossOfFunction(changeEffects);
				boolean islof = lof.isLof();
				Assert.assertEquals(true, islof);
			}
		}

	}

	/**
	 * Check that Core Splice Donor acceptors are considered LOF
	 * @param tr
	 * @param ex
	 * @param exonNum
	 */
	void checkSpliceDonor(Transcript tr, Exon ex, int exonNum) {
		int step = tr.isStrandPlus() ? 1 : -1;
		int maxRank = tr.numChilds();
		int intronNum = exonNum; // We care about the intron before

		if (ex.getRank() < maxRank) {
			// Position
			int posDonor = tr.isStrandPlus() ? ex.getEnd() : ex.getStart();
			posDonor += step;

			// Splice site size
			int maxSize = Math.min(tr.intronSize(intronNum), SpliceSite.CORE_SPLICE_SITE_SIZE);
			if (debug) Gpr.debug("Intron size: " + tr.intronSize(intronNum));
			if (maxSize <= 0) throw new RuntimeEOFException("Max splice size is " + maxSize);

			//---
			// For all position on splice site donor positions, make sure it is LOF
			//---
			for (int pos = posDonor, i = 0; i < maxSize; i++, pos += step) {
				SeqChange seqChange = new SeqChange(tr.getChromosome(), pos, pos, ""); // Create a seqChange
				Marker marker = findMarker(seqChange, EffectType.SPLICE_SITE_DONOR, null, ex);
				LinkedList<ChangeEffect> changeEffects = changeEffects(seqChange, EffectType.SPLICE_SITE_DONOR, marker); // Create a SPLICE_DONOR effect
				if (debug) Gpr.debug("SeqChange:" + seqChange);

				// Create a LOF object and analyze the effect
				LossOfFunction lof = new LossOfFunction(changeEffects);
				boolean islof = lof.isLof();
				Assert.assertEquals(true, islof);
			}
		}

	}

	/**
	 * Find a marker that intersects seqChange
	 * @return
	 */
	Marker findMarker(SeqChange seqChange, EffectType effectType, Transcript tr, Exon exon) {
		Markers markers = config.getSnpEffectPredictor().intersects(seqChange);
		for (Marker m : markers) {
			Exon mex = (Exon) m.findParent(Exon.class);
			Transcript mtr = (Transcript) m.findParent(Transcript.class);

			if ((m.getType() == effectType) && (mex != null) && (mtr != null)) {
				if (exon != null) {
					// Exon filter?
					if (mex.getId().equals(exon.getId())) return m;
				} else if (tr != null) {
					// Transcript filter?
					if (mtr.getId().equals(tr.getId())) return m;
				} else return m; // No exon reference? => just return this
			}
		}

		throw new RuntimeEOFException("Cannot find '" + effectType + "' " + (exon != null ? "for exon " + exon.getId() : "") + ", seqChange: " + seqChange);
	}

	public void test_01() {
		// Load database
		String genomeVer = "testHg3766Chr1";
		Gpr.debug("Loading database '" + genomeVer + "'");
		config = new Config(genomeVer, Config.DEFAULT_CONFIG_FILE);
		config.loadSnpEffectPredictor();
		Gpr.debug("Building forest");
		config.getSnpEffectPredictor().buildForest();

		// For each gene, transcript, check that NMD works
		Gpr.debug("Testing");
		for (Gene gene : config.getGenome().getGenes()) {
			System.err.println("LOF test\tGene ID:" + gene.getId());
			for (Transcript tr : gene) {
				//				if (tr.getId().equals("ENST00000327057")) {
				if (debug) System.err.println(tr);
				checkLof(tr);
				//				}
			}
		}
	}
}
