package ca.mcgill.mcb.pcingola.snpEffect.testCases;

import java.util.LinkedList;

import junit.framework.TestCase;
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

	Config config;

	public TestCasesNmd() {
		super();
	}

	public void test_01() {
		// Load database
		String genomeVer = "testHg3766Chr1";
		Gpr.debug("Loading database '" + genomeVer + "'");
		Config config = new Config(genomeVer, Config.DEFAULT_CONFIG_FILE);
		config.loadSnpEffectPredictor();

		// For each gene, transcript, exon, position
		for (Gene gene : config.getGenome().getGenes()) {
			for (Transcript tr : gene) {

				//				if (tr.getId().equals("ENST00000367789")) {
				//for (Exon ex : tr) {
				// for (int pos = ex.getStart(); pos < ex.getEnd(); pos++) {

				int pos = tr.getStart();
				// Gpr.debug(gene.getId() + "\t" + tr.getId() + "\t" + ex.getId() + "\t" + pos);
				Gpr.debug(gene.getId() + "\t" + tr.getId() + "\t" + pos);

				// Create a seqChange
				SeqChange seqChange = new SeqChange(tr.getChromosome(), pos, pos, "");

				// Create a STOP_GAIN effect
				ChangeEffect changeEffect = new ChangeEffect(seqChange);
				//changeEffect.set(ex, EffectType.STOP_GAINED, "");
				changeEffect.set(tr, EffectType.STOP_GAINED, "");
				LinkedList<ChangeEffect> changeEffects = new LinkedList<ChangeEffect>();
				changeEffects.add(changeEffect);

				// Create a LOF object and analyze the effect
				LossOfFunction lof = new LossOfFunction(changeEffects);
				//				lof.isNmd();
				lof.lastNmdPos(tr);
				// }
				// }
				//				}
			}
		}
	}
}
