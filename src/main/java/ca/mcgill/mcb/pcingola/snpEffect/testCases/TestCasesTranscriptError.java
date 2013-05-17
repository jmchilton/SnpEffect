package ca.mcgill.mcb.pcingola.snpEffect.testCases;

import java.util.List;

import junit.framework.Assert;
import junit.framework.TestCase;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectImpact;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.WarningType;
import ca.mcgill.mcb.pcingola.snpEffect.commandLine.SnpEffCmdEff;
import ca.mcgill.mcb.pcingola.vcf.VcfEffect;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;

/**
 * Test case where VCF entries hit a transcript that has errors
 * 
 * @author pcingola
 */
public class TestCasesTranscriptError extends TestCase {

	public TestCasesTranscriptError() {
		super();
	}

	public void test_01() {
		String args[] = { "testHg3763Chr20", "./tests/short_codon_bug.vcf" };
		transcriptError(args, WarningType.WARNING_TRANSCRIPT_INCOMPLETE);
	}

	/**
	 * Run a predictor and check if the expected warnings appear
	 * @param args
	 * @param warningType
	 */
	void transcriptError(String args[], WarningType warningType) {
		SnpEffCmdEff snpEffCmdEff = new SnpEffCmdEff();
		snpEffCmdEff.parseArgs(args);
		snpEffCmdEff.setVerbose(true);
		List<VcfEntry> vcfEntries = snpEffCmdEff.run(true);

		boolean hasWarning = false;
		for (VcfEntry ve : vcfEntries) {
			System.out.println(ve);
			for (VcfEffect veff : ve.parseEffects()) {
				EffectImpact imp = veff.getImpact();
				System.out.println("\t" + imp + "\t" + veff);

				// Chek if the warning type we expect is there
				if (veff.getErrorsOrWarning() != null) hasWarning |= veff.getErrorsOrWarning().indexOf(warningType.toString()) >= 0;
			}
		}

		Assert.assertEquals(true, hasWarning);
	}
}
