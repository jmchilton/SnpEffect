package ca.mcgill.mcb.pcingola.snpEffect.commandLine;

import java.util.HashSet;

import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.snpEffect.LossOfFunction;
import ca.mcgill.mcb.pcingola.snpEffect.LossOfFunctionEntry;
import ca.mcgill.mcb.pcingola.snpEffect.NonsenseMediatedDecayEntry;
import ca.mcgill.mcb.pcingola.stats.CountByType;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;
import ca.mcgill.mcb.pcingola.vcf.VcfEffect;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;

/**
 * Command line: Test
 * 
 * Note: Used for testing weird stuff
 * 
 * @author pcingola
 */
public class SnpEffCmdTest extends SnpEff {

	public static final String VARIANTS_IN_GENES = "_VARAINTS_IN_GENES";
	public static final String VARIANTS = "_VARAINTS";
	public static final int SHOW_EVERY = 10000;

	String genesFile;
	String vcfFile;
	CountByType countByEff;
	HashSet<String> genes = new HashSet<String>();

	public SnpEffCmdTest() {
		super();
		genes = new HashSet<String>();
		countByEff = new CountByType();
	}

	/**
	 * Analyze vcf entry
	 * @param ve
	 */
	void analyze(VcfEntry ve) {
		boolean inGenes = false;

		//---
		// Parse Effect
		//---
		for (VcfEffect veff : ve.parseEffects()) {

			// Do not process is there are errors or warnings
			if (veff.getErrorsOrWarning() != null) continue;

			// No gene info? Nothing to do
			String gene = veff.getGene();
			if (gene == null || gene.isEmpty()) continue;

			// Gene Info does not match? Nothing to do
			if (!genes.contains(gene)) continue;

			inGenes = true;

			// Count by effect
			countByEff.inc(veff.getEffect().toString());
			if (veff.getEffectDetails() != null && !veff.getEffectDetails().isEmpty()) countByEff.inc(veff.getEffect() + "[" + veff.getEffectDetails() + "]");
		}

		//---
		// Parse LOF 
		//---
		for (LossOfFunctionEntry lof : ve.parseLof()) {
			// No gene info? Nothing to do
			String gene = lof.getGeneName();
			if (gene == null || gene.isEmpty()) continue;

			// Gene Info does not match? Nothing to do
			if (!genes.contains(gene)) continue;

			inGenes = true;
			countByEff.inc(LossOfFunction.VCF_INFO_LOF_NAME);
		}

		//---
		// Parse NMD 
		//---
		for (NonsenseMediatedDecayEntry nmd : ve.parseNmd()) {
			// No gene info? Nothing to do
			String gene = nmd.getGeneName();
			if (gene == null || gene.isEmpty()) continue;

			// Gene Info does not match? Nothing to do
			if (!genes.contains(gene)) continue;

			inGenes = true;
			countByEff.inc(LossOfFunction.VCF_INFO_NMD_NAME);
		}

		// Count if it is in genes
		if (inGenes) countByEff.inc(VARIANTS_IN_GENES);

		// Count total number of variants
		countByEff.inc(VARIANTS);
	}

	/**
	 * Parse command line arguments
	 */
	@Override
	public void parseArgs(String[] args) {
		if (args.length != 2) usage("Missing arguments");
		vcfFile = args[0];
		genesFile = args[1];
	}

	/**
	 * Run command
	 */
	@Override
	public boolean run() {
		//---
		// Load genes
		//---
		if (verbose) Timer.showStdErr("Loading genes from '" + genesFile + "'");
		for (String gene : Gpr.readFile(genesFile).split("\n"))
			genes.add(gene.trim());
		if (verbose) Timer.showStdErr("Done. Genes added : " + genes.size());

		//---
		// Count effect in VCF
		//---
		if (verbose) Timer.showStdErr("Counting effect on input: " + vcfFile);
		VcfFileIterator vcf = new VcfFileIterator(vcfFile);
		int i = 1;
		for (VcfEntry ve : vcf) {
			analyze(ve);

			if (verbose) Gpr.showMark(i++, SHOW_EVERY);
		}

		System.out.println("GENES\t" + genes.size());
		System.out.println(countByEff);
		return true;
	}

	/**
	 * Show usage and exit
	 */
	@Override
	public void usage(String message) {
		if (message != null) System.err.println("Error: " + message + "\n");
		System.err.println("snpEff version " + SnpEff.VERSION);
		System.err.println("Usage: snpEff test file.vcf genes.txt");
		System.exit(-1);
	}

}
