package ca.mcgill.mcb.pcingola.snpEffect.testCases;

import java.util.List;

import junit.framework.TestCase;

import org.junit.Assert;

import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectImpact;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;
import ca.mcgill.mcb.pcingola.snpEffect.commandLine.SnpEffCmdEff;
import ca.mcgill.mcb.pcingola.vcf.VcfEffect;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;

/**
 * Test Motif databases
 * 
 * @author pcingola
 */
public class TestCasesMotif extends TestCase {

	public static boolean debug = false;
	public static boolean verbose = true;
	public static int SHOW_EVERY = 10;

	public TestCasesMotif() {
		super();
	}

	void checkMotif(String genomeVer, String vcfFile, String effectDetails, EffectImpact impact) {
		String args[] = { "-motif", "-ud", "0", genomeVer, vcfFile };
		SnpEffCmdEff cmd = new SnpEffCmdEff();
		cmd.setVerbose(true);
		cmd.parseArgs(args);

		// Run
		List<VcfEntry> vcfEntries = cmd.run(true);

		// Check results
		int numNextProt = 0;
		for (VcfEntry ve : vcfEntries) {
			for (VcfEffect veff : ve.parseEffects()) {
				System.out.println("\t" + veff);
				if ((veff.getEffect() == EffectType.MOTIF) // Is it motif?
						&& effectDetails.equals(veff.getEffectDetails()) // Are details OK? 					
						&& (impact == veff.getImpact())) // Is impact OK?
					numNextProt++;
			}
		}

		Assert.assertEquals(1, numNextProt);
	}

	public void test_01() {
		checkMotif("testHg3770Chr22", "tests/test_motif_01.vcf", "MA0099.2:AP1", EffectImpact.LOW);
	}

	public void test_02() {
		checkMotif("testHg3770Chr22", "tests/test_motif_02.vcf", "MA0099.2:AP1", EffectImpact.MODIFIER);
	}

	public void test_03() {
		checkMotif("testHg3770Chr22", "tests/test_motif_03.vcf", "MA0099.2:AP1", EffectImpact.LOW);
	}

}
