package ca.mcgill.mcb.pcingola.snpEffect.testCases;

import junit.framework.TestCase;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * Test Loss of Function prediction
 * 
 * @author pcingola
 */
public class TestCasesLof extends TestCase {

	public static boolean debug = false;

	Config config;

	public TestCasesLof() {
		super();
	}

	/**
	 * Check that LOF works for a given transcript
	 * @param gene
	 * @param tr
	 */
	void checkLof(Gene gene, Transcript tr) {
		checkLofExonDeleted(tr);
		checkLofFrameShift(tr);
		checkLofSplice(tr);
		checkLofStartLost(tr);
	}

	void checkLofExonDeleted(Transcript tr) {
		throw new RuntimeException("Unimplemented!");
	}

	void checkLofFrameShift(Transcript tr) {
		throw new RuntimeException("Unimplemented!");
	}

	void checkLofSplice(Transcript tr) {
		throw new RuntimeException("Unimplemented!");
	}

	void checkLofStartLost(Transcript tr) {
		throw new RuntimeException("Unimplemented!");
	}

	public void test_01() {
		// Load database
		String genomeVer = "testHg3766Chr1";
		Gpr.debug("Loading database '" + genomeVer + "'");
		Config config = new Config(genomeVer, Config.DEFAULT_CONFIG_FILE);
		config.loadSnpEffectPredictor();

		// For each gene, transcript, check that NMD works
		for (Gene gene : config.getGenome().getGenes()) {
			System.err.println("NMD test\tGene ID:" + gene.getId());
			for (Transcript tr : gene) {
				//				if (tr.getId().equals("ENST00000486084")) {
				if (debug) System.err.println(tr);
				checkLof(gene, tr);
				//				}
			}
		}
	}
}
