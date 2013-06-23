package ca.mcgill.mcb.pcingola.snpEffect.commandLine;

import java.util.HashSet;

import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.LossOfFunction;
import ca.mcgill.mcb.pcingola.snpEffect.LossOfFunctionEntry;
import ca.mcgill.mcb.pcingola.snpEffect.NonsenseMediatedDecayEntry;
import ca.mcgill.mcb.pcingola.snpEffect.SnpEffectPredictor;
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
	public static final double MIN_PERCENT_TRANSCRIPTS_AFFECTED = 0.0;

	public static final int SHOW_EVERY = 10000;

	boolean onlyProteinCodingTranscripts = false; // Use only protein coding transcripts
	boolean useClosestGene = true;
	boolean doNotUseAF50 = false;

	SnpEffectPredictor snpEffectPredictor;
	String genesFile;
	String vcfFile;
	CountByType count; // Counter
	CountByType countByEffect; // Count raw number of effects
	CountByType countByVariant; // Count each effect once per variant
	CountByType countByEffByGene; // Count number of effects for each gene
	CountByType countByGene;// Count each effect once per gene
	HashSet<String> genes; // Only select effect in these genes

	public SnpEffCmdTest() {
		super();
		countByEffByGene = new CountByType();
		countByVariant = new CountByType();
		countByGene = new CountByType();
		countByEffect = new CountByType();
		count = new CountByType();
	}

	/**
	 * Analyze vcf entry
	 * @param ve
	 */
	void analyze(VcfEntry ve) {
		boolean inGenes = false;
		HashSet<String> effectsByVariant = new HashSet<String>();
		HashSet<String> effectsByGene = new HashSet<String>();

		//		// We ignore AF > 0.5 because 
		//		// Most eQtl algorithms are based on minor allele frequencies. So when 
		//		// eQtl is calculated, REF and ALT are swapped.
		//		// This means that the EFF field is actually not describing the effect of the eQTL and we should filter it out
		//		if (ve.getInfoFloat("AF") > 0.5) return;

		String geneClosest = useClosestGene ? findClosestGene(ve) : "";

		//---
		// Parse Effect
		//---
		for (VcfEffect veff : ve.parseEffects()) {

			// Do not process is there are errors or warnings
			if (veff.getErrorsOrWarning() != null) {
				count.inc("ERRORS_OR_WARNINGS");
				continue;
			}

			String gene = veff.getGene();
			if (genes != null) {
				// No gene info? Nothing to do
				if (gene == null || gene.isEmpty()) {
					count.inc("NO_GENE");
					continue;
				}

				// Gene Info does not match? Nothing to do
				if (!genes.contains(gene)) {
					count.inc("NO_GENE_SET");
					continue;
				}

				// Not a protein coding transcript? Skip
				if (onlyProteinCodingTranscripts && ((veff.getBioType() == null) || !veff.getBioType().equals("protein_coding"))) {
					count.inc(BIOTYPE_SKIPPED + "_" + veff.getBioType());
					continue;
				}

				inGenes = true;
			} else if (gene == null || gene.isEmpty()) {
				if (!geneClosest.isEmpty()) count.inc("GENE_CLOSEST");
				gene = geneClosest; // Use closest gene 
			}

			// Count by effect
			String key = veff.getEffect().toString();
			if (veff.getEffectDetails() != null && !veff.getEffectDetails().isEmpty()) key += "[" + veff.getEffectDetails() + "]";
			effectsByVariant.add(key);
			effectsByGene.add(gene + "\t" + key);
			countByEffect.inc(key);
		}

		//---
		// Parse LOF 
		//---
		for (LossOfFunctionEntry lof : ve.parseLof()) {
			// No gene info? Nothing to do
			String gene = lof.getGeneName();
			if (gene == null || gene.isEmpty()) continue;

			// Gene Info does not match? Nothing to do
			if ((genes != null) && !genes.contains(gene)) continue;

			inGenes = true;
			if (lof.getPercentOfTranscriptsAffected() >= MIN_PERCENT_TRANSCRIPTS_AFFECTED) {
				effectsByGene.add(gene + "\t" + LossOfFunction.VCF_INFO_LOF_NAME);
				effectsByVariant.add(LossOfFunction.VCF_INFO_LOF_NAME);
				countByEffect.inc(LossOfFunction.VCF_INFO_LOF_NAME);
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
			if ((genes != null) && !genes.contains(gene)) continue;

			inGenes = true;
			if (nmd.getPercentOfTranscriptsAffected() >= MIN_PERCENT_TRANSCRIPTS_AFFECTED) {
				effectsByGene.add(gene + "\t" + LossOfFunction.VCF_INFO_NMD_NAME);
				effectsByVariant.add(LossOfFunction.VCF_INFO_NMD_NAME);
				countByEffect.inc(LossOfFunction.VCF_INFO_NMD_NAME);
			}
		}

		// ACAT & NCCAT Scores (they appear once per variant)
		String acat = ve.getInfo(SnpEffCmdAcat.ACAT);
		if (acat != null) {
			String acatFields[] = acat.split(",");
			for (String af : acatFields) {
				String afs[] = af.split(":");
				effectsByVariant.add("_" + SnpEffCmdAcat.ACAT + "_" + afs[2]);
			}
		}

		// NCCAT is just once per variant
		String nccat = ve.getInfo(SnpEffCmdAcat.NCCAT);
		if (nccat != null) countByVariant.inc("_" + SnpEffCmdAcat.NCCAT + "_" + nccat);

		// Count once per variant
		for (String eff : effectsByVariant)
			countByVariant.inc(eff);

		// Count effects by gene
		for (String eff : effectsByGene)
			countByEffByGene.inc(eff);

		// Count total number of variants
		count.inc(VARIANTS);
		if (inGenes) count.inc(VARIANTS_IN_GENES); // Count if it is in genes
	}

	/**
	 * Find closes gene name
	 * @param queryMarker
	 * @return
	 */
	String findClosestGene(Marker queryMarker) {
		Gene gene = snpEffectPredictor.queryClosestGene(queryMarker);
		return gene != null ? gene.getId() : "";
	}

	/**
	 * Parse command line arguments
	 */
	@Override
	public void parseArgs(String[] args) {
		if (args.length < 2) usage(null);

		int idx = 0;
		genomeVer = args[idx++];
		vcfFile = args[idx++];
		if (args.length > idx) genesFile = args[idx];
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
		// Load database, build tree
		//---
		if (verbose) Timer.showStdErr("Reading configuration...");
		config = new Config(genomeVer, configFile); // Read configuration
		if (verbose) Timer.showStdErr("done");

		if (verbose) Timer.showStdErr("Loading predictor...");
		config.loadSnpEffectPredictor();
		if (verbose) Timer.showStdErr("done");

		if (verbose) Timer.showStdErr("Building interval forest...");
		snpEffectPredictor = config.getSnpEffectPredictor();
		snpEffectPredictor.buildForest();
		if (verbose) Timer.showStdErr("done");

		//---
		// Load genes
		//---
		if (genesFile != null) {
			genes = new HashSet<String>();
			if (verbose) Timer.showStdErr("Loading genes from '" + genesFile + "'");
			for (String gene : Gpr.readFile(genesFile).split("\n"))
				genes.add(gene.trim());
			if (verbose) Timer.showStdErr("Done. Genes added : " + genes.size());
		}

		//---
		// Process input file
		//---
		if (verbose) Timer.showStdErr("Counting effect on input: " + vcfFile);
		VcfFileIterator vcf = new VcfFileIterator(vcfFile);
		int i = 1;
		for (VcfEntry ve : vcf) {
			analyze(ve);

			if (verbose) Gpr.showMark(i++, SHOW_EVERY);
		}

		//---
		// Calculate 'once per gene' counters
		//---
		for (String key : countByEffByGene.keySet()) {
			if (countByEffByGene.get(key) > 0) { // This should always be true
				String keySplit[] = key.split("\t");

				if (keySplit.length > 1) {
					String eff = keySplit[1];
					countByGene.inc(eff);
				}
			} else throw new RuntimeException("This should never happen!");
		}

		//---
		// Show output
		//---
		System.out.println("# General Numbers");
		if (genes != null) System.out.println("GENES\t" + genes.size());
		print("", count);
		System.out.println("#");
		System.out.println("# Number of effects (raw counts)");
		print("count_effect", countByEffect);
		System.out.println("#");
		System.out.println("# Number of effects per gene (number of effects for each gene)");
		print("count_effect_by_gene", countByEffByGene);
		System.out.println("#");
		System.out.println("# Number of effects per variant (i.e. each effect is counted only once per variant)");
		print("count_by_variant", countByVariant);
		System.out.println("#");
		System.out.println("# Number of genes for each effect (i.e. each effect is counted only once per gene)");
		print("count_by_gene", countByGene);
		return true;
	}

	/**
	 * Show usage and exit
	 */
	@Override
	public void usage(String message) {
		if (message != null) System.err.println("Error: " + message + "\n");
		System.err.println("snpEff version " + SnpEff.VERSION);
		System.err.println("Usage: snpEff test genomeVer file.vcf [genes.txt]");
		System.exit(-1);
	}

}
