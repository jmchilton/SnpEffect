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
	public static final String BIOTYPE_SKIPPED = "_BIOTYPE_SKIPPED";
	public static final double MIN_PERCENT_TRANSCRIPTS_AFFECTED = 0.5;

	public static final int SHOW_EVERY = 10000;

	String genesFile;
	String vcfFile;
	CountByType countByGene;
	CountByType countByVariant;
	HashSet<String> genes = new HashSet<String>();

	public SnpEffCmdTest() {
		super();
		genes = new HashSet<String>();
		countByGene = new CountByType();
		countByVariant = new CountByType();
	}

	/**
	 * Analyze vcf entry
	 * @param ve
	 */
	void analyze(VcfEntry ve) {
		boolean inGenes = false;
		HashSet<String> effectsByVariant = new HashSet<String>();

		// We ignore AF > 0.5 because 
		// Most eQtl algorithms are based on minor allele frequencies. So when 
		// eQtl is calculated, REF and ALT are swapped.
		// This means that the EFF field is actually not describing the effect of the eQTL and we should filter it out
		if (ve.getInfoFloat("AF") > 0.5) return;

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

			// Not a protrein coding transcript? Skip
			if ((veff.getBioType() == null) || !veff.getBioType().equals("protein_coding")) {
				countByGene.inc(BIOTYPE_SKIPPED + "_" + veff.getBioType());
				continue;
			}

			inGenes = true;

			// Count by effect
			if (veff.getEffectDetails() != null && !veff.getEffectDetails().isEmpty()) {
				countByGene.inc(gene + "\t" + veff.getEffect() + "[" + veff.getEffectDetails() + "]");
				effectsByVariant.add(veff.getEffect() + "[" + veff.getEffectDetails() + "]");
			} else {
				countByGene.inc(gene + "\t" + veff.getEffect().toString());
				effectsByVariant.add(veff.getEffect().toString());
			}
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
			if (lof.getPercentOfTranscriptsAffected() >= MIN_PERCENT_TRANSCRIPTS_AFFECTED) {
				countByGene.inc(gene + "\t" + LossOfFunction.VCF_INFO_LOF_NAME);
				effectsByVariant.add(LossOfFunction.VCF_INFO_LOF_NAME);
			}
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
			if (nmd.getPercentOfTranscriptsAffected() >= MIN_PERCENT_TRANSCRIPTS_AFFECTED) {
				countByGene.inc(gene + "\t" + LossOfFunction.VCF_INFO_NMD_NAME);
				effectsByVariant.add(LossOfFunction.VCF_INFO_NMD_NAME);
			}
		}

		// Count once per variant
		for (String eff : effectsByVariant)
			countByVariant.inc(eff);

		// Count if it is in genes
		if (inGenes) countByGene.inc(VARIANTS_IN_GENES);

		// Count total number of variants
		countByGene.inc(VARIANTS);
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

	void print(String label, CountByType countByType) {
		System.out.println(label + "\teff\tcount");
		for (String type : countByType.keysSorted())
			System.out.println(label + "\t" + type + "\t" + countByType.get(type));

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

		//---
		// Increment all counters. Only once per gene
		//---
		CountByType countByEff = new CountByType(); // Count each effect only once per gene
		for (String key : countByGene.keySet()) {
			if (countByGene.get(key) > 0) { // This should always be true
				String keySplit[] = key.split("\t");

				if (keySplit.length > 1) {
					String eff = keySplit[1];
					countByEff.inc(eff);
				}
			} else throw new RuntimeException("This should never happen!");
		}

		System.out.println("# General Numbers");
		System.out.println("GENES\t" + genes.size());
		System.out.println("VARIANTS\t" + countByGene.get(VARIANTS));
		System.out.println("VARIANTS_IN_GENES\t" + countByGene.get(VARIANTS_IN_GENES));
		System.out.println("#");
		System.out.println("# Number of effects per gene (number of effects for each gene)");
		print("count_effect_by_gene", countByGene);
		System.out.println("#");
		System.out.println("# Number of effects per variant (i.e. each effect is counted only once per variant)");
		print("count_by_variant", countByVariant);
		System.out.println("#");
		System.out.println("# Number of genes for each effect (i.e. each effect is counted only once per gene)");
		print("count_by_gene", countByEff);
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
