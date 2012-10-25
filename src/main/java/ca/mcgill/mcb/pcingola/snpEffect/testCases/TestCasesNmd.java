package ca.mcgill.mcb.pcingola.snpEffect.testCases;

import java.util.HashSet;
import java.util.LinkedList;

import junit.framework.TestCase;
import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.LossOfFunction;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * Test Nonsense mediated decay prediction
 * 
 * @author pcingola
 */
public class TestCasesNmd extends TestCase {

	public static boolean debug = false;

	Config config;

	public TestCasesNmd() {
		super();
	}

	/**
	 * Check that NMD works for a given transcript
	 * @param gene
	 * @param tr
	 */
	void checkNmd(Gene gene, Transcript tr) {
		System.err.print("\tTranscript " + tr.getId() + " " + (tr.isStrandPlus() ? '+' : '-') + " :");
		int pos = 0;
		boolean isNmd[] = new boolean[tr.cds().length()];
		Exon latestCodingExon = null;
		HashSet<Exon> codingExons = new HashSet<Exon>();

		for (Exon exon : tr.sortedStrand()) {
			for (int expos = exon.getStart(); expos < exon.getEnd(); expos++) {
				// Not in UTR? => Test
				if (!tr.isUtr(expos)) {
					latestCodingExon = exon;
					codingExons.add(exon);

					// Create a seqChange
					SeqChange seqChange = new SeqChange(tr.getChromosome(), expos, expos, "");

					// Create a STOP_GAIN effect
					ChangeEffect changeEffect = new ChangeEffect(seqChange);
					changeEffect.set(exon, EffectType.STOP_GAINED, "");
					LinkedList<ChangeEffect> changeEffects = new LinkedList<ChangeEffect>();
					changeEffects.add(changeEffect);

					// Create a LOF object and analyze the effect
					LossOfFunction lof = new LossOfFunction(changeEffects);
					int lastNmdPos = lof.lastNmdPos(tr);
					isNmd[pos] = lof.isNmd();

					// Show info
					if (debug) Gpr.debug(gene.getId() + "\t" + tr.getId() + "\t" + exon.getId() //
							+ "\n\tisNmd[" + pos + "]      : " + isNmd[pos] // 
							+ "\n\tlastNmdPos   : " + lastNmdPos // 
							+ "\n\tSeqChange    : " + seqChange //
							+ "\n\tChangeEffect : " + changeEffect //
					);

					System.err.print(isNmd[pos] ? '+' : '.');
					pos++;
				} else System.err.print('U');
			}
			System.err.print(" | ");
		}
		System.err.println("");

		//---
		// Check that NMP prediction is 'correct'
		//---
		// We need a spilce event in the coding part 
		if (codingExons.size() > 1) {
			if (latestCodingExon == null) throw new RuntimeException("Cannot find latest coding exon!");
		}
	}

	public void test_01() {
		// Load database
		String genomeVer = "testHg3766Chr1";
		Gpr.debug("Loading database '" + genomeVer + "'");
		Config config = new Config(genomeVer, Config.DEFAULT_CONFIG_FILE);
		config.loadSnpEffectPredictor();

		// For each gene, transcript, check that NMD works
		for (Gene gene : config.getGenome().getGenes()) {
			//System.err.println("Gene ID:" + gene.getId());
			for (Transcript tr : gene) {
				if (tr.getId().equals("ENST00000422560")) checkNmd(gene, tr);
			}
		}
	}
}
