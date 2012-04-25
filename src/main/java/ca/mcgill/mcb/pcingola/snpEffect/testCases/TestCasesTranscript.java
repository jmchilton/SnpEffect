package ca.mcgill.mcb.pcingola.snpEffect.testCases;

import java.util.Random;

import junit.framework.Assert;
import junit.framework.TestCase;
import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Genome;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.SnpEffectPredictor;
import ca.mcgill.mcb.pcingola.snpEffect.factory.SnpEffPredictorFactoryRand;

/**
 * Test random SNP changes 
 * 
 * @author pcingola
 */
public class TestCasesTranscript extends TestCase {

	boolean debug = false;
	Random rand;
	Config config;
	Genome genome;
	Chromosome chromosome;
	Gene gene;
	Transcript transcript;
	SnpEffectPredictor snpEffectPredictor;
	String chromoSequence = "";
	char chromoBases[];

	public TestCasesTranscript() {
		super();
		init();
	}

	void init() {
		initRand();
		initSnpEffPredictor();
	}

	void initRand() {
		rand = new Random(20120131);
	}

	void initSnpEffPredictor() {
		// Create a config and force out snpPredictor for hg37 chromosome Y
		config = new Config("testCase", Config.DEFAULT_CONFIG_FILE);

		// Create factory
		int maxGeneLen = 1000;
		int maxTranscripts = 1;
		int maxExons = 5;
		SnpEffPredictorFactoryRand sepf = new SnpEffPredictorFactoryRand(config, 1, rand, maxGeneLen, maxTranscripts, maxExons);

		// Create predictor
		snpEffectPredictor = sepf.create();
		config.setSnpEffectPredictor(snpEffectPredictor);

		// Chromosome sequence
		chromoSequence = sepf.getChromoSequence();
		chromoBases = chromoSequence.toCharArray();

		// No upstream or downstream
		config.getSnpEffectPredictor().setUpDownStreamLength(0);

		// Build forest
		config.getSnpEffectPredictor().buildForest();

		chromosome = sepf.getChromo();
		genome = config.getGenome();
		gene = genome.getGenes().iterator().next();
		transcript = gene.iterator().next();
	}

	public void test_CdsPos() {
		int N = 1000;

		// Test N times: 
		//		- Create a random gene transcript, exons
		// 		- Cal
		for( int iter = 0; iter < N; iter++ ) {
			initSnpEffPredictor();
			if( debug ) System.out.println("Test CDS pos iteration: " + iter + "\n" + transcript);
			else System.out.println("Test CDS pos iteration: " + iter + "\t" + (transcript.getStrand() >= 0 ? "+" : "-") + "\t" + transcript.cds());

			int cdsBaseNum = 0;
			int cds2pos[] = transcript.cdsBaseNumber2ChrPos();

			// For each exon...
			for( Exon exon : transcript.sortedStrand() ) {
				// Iterate on each base and compare CDS positon with calculated one
				int min = transcript.isStrandPlus() ? exon.getStart() : exon.getEnd();
				int step = transcript.isStrandPlus() ? 1 : -1;
				for( int pos = min; exon.intersects(pos); pos += step, cdsBaseNum++ ) {
					int cdsBaseNumCalc = transcript.cdsBaseNumber(pos, true);

					// Is it OK?
					Assert.assertEquals(cdsBaseNum, cdsBaseNumCalc);
					Assert.assertEquals(pos, cds2pos[cdsBaseNum]);
				}
			}
		}
	}
}
