package ca.mcgill.mcb.pcingola.snpEffect.commandLine;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import ca.mcgill.mcb.pcingola.collections.AutoHashMap;
import ca.mcgill.mcb.pcingola.fileIterator.BedFileIterator;
import ca.mcgill.mcb.pcingola.fileIterator.LineFileIterator;
import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.geneSets.GeneSets;
import ca.mcgill.mcb.pcingola.geneSets.GeneSetsRanked;
import ca.mcgill.mcb.pcingola.geneSets.algorithm.EnrichmentAlgorithm;
import ca.mcgill.mcb.pcingola.geneSets.algorithm.EnrichmentAlgorithm.EnrichmentAlgorithmType;
import ca.mcgill.mcb.pcingola.geneSets.algorithm.EnrichmentAlgorithmGreedyVariableSize;
import ca.mcgill.mcb.pcingola.geneSets.algorithm.FisherPValueAlgorithm;
import ca.mcgill.mcb.pcingola.geneSets.algorithm.FisherPValueGreedyAlgorithm;
import ca.mcgill.mcb.pcingola.geneSets.algorithm.LeadingEdgeFractionAlgorithm;
import ca.mcgill.mcb.pcingola.geneSets.algorithm.RankSumPValueAlgorithm;
import ca.mcgill.mcb.pcingola.geneSets.algorithm.RankSumPValueGreedyAlgorithm;
import ca.mcgill.mcb.pcingola.gsa.ChrPosPvalueList;
import ca.mcgill.mcb.pcingola.gsa.PvalueList;
import ca.mcgill.mcb.pcingola.gsa.PvalueList.PvalueSummary;
import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Genome;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Markers;
import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.SnpEffectPredictor;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;

/**
 * Command line: Gene-Sets Analysis
 * 
 * Perform gene set analysys
 * 
 * @author pcingola
 */
public class SnpEffCmdGsa extends SnpEff {

	public enum CorrectionMethod {
		NONE
	}

	public static int READ_INPUT_SHOW_EVERY = 1000;;

	InputFormat inputFormat = InputFormat.VCF;
	boolean useClosestGene = false; // Map to 'any' closest gene?
	boolean useGeneId = false; // Use geneId instead of geneName
	int upDownStreamLength = SnpEffectPredictor.DEFAULT_UP_DOWN_LENGTH;
	int minGeneSetSize = 0;
	int maxGeneSetSize = Integer.MAX_VALUE;
	int numberofGeneSetsToSelect = Integer.MAX_VALUE; // TODO : Add a command line option to limit this
	int initGeneSetSize = 100;
	int randIterations = 0;
	double maxPvalue = 1.0;
	double interestingPerc = 0.95;
	String inputFile = "";
	String infoName = "";
	String msigdb = "";
	String genePvalueFile = "";
	PvalueSummary pvalueSummary = PvalueSummary.MIN;
	CorrectionMethod correctionMethod = CorrectionMethod.NONE;
	SnpEffectPredictor snpEffectPredictor;
	Genome genome;
	GeneSets geneSets;
	ChrPosPvalueList chrPosPvalueList; // List of <chr,pos,pvalue>
	AutoHashMap<String, PvalueList> genePvalues; // A map of geneId -> List[pValues]
	HashMap<String, Double> genePvalue; // A <gene, pValue> map
	EnrichmentAlgorithmType enrichmentAlgorithmType = EnrichmentAlgorithmType.RANKSUM_GREEDY;

	public SnpEffCmdGsa() {
		super();
	}

	/**
	 * Correct p-values (e.g. using covariates)
	 */
	void correctPvalues() {
		switch (correctionMethod) {
		case NONE:
			// Nothing to do
			break;

		default:
			throw new RuntimeException("Unimplemented p-value correction method '" + correctionMethod + "'");
		}
	}

	/**
	 * Create interesting genes
	 */
	void createInterestingGenes() {
		// Get 
		PvalueList pvalues = new PvalueList();
		for (double pval : genePvalue.values())
			pvalues.add(pval);

		// Get p-value threshold
		double pThreshold = pvalues.quantile(1.0 - interestingPerc);

		// Mark all p-values lower than that as 'interesting'
		int count = 0, countAdded = 0;
		geneSets.setDoNotAddIfNotInGeneSet(true);
		for (String geneId : genePvalue.keySet())
			if (genePvalue.get(geneId) <= pThreshold) {
				if (geneSets.addInteresting(geneId)) countAdded++; // Count added genes
				count++;
			}

		// Show info
		if (verbose) {
			double realPerc = (100.0 * count) / genePvalue.size();
			double realPercAdded = (100.0 * countAdded) / genePvalue.size();
			Timer.showStdErr(String.format("P-value threshold:"//
					+ "\n\tQuantile                 : %.2f%%"//
					+ "\n\tThreshold                : %f"//
					+ "\n\tInteresting genes        : %d  (%.2f%%)" //
					+ "\n\tInteresting genes  added : %d  (%.2f%%)" //
			, 100.0 * interestingPerc, pThreshold, count, realPerc, countAdded, realPercAdded));
		}
	}

	/**
	 * Perform enrichment analysis
	 */
	void enrichmentAnalysis() {
		// Initialize gene set values
		for (String geneId : genePvalue.keySet())
			geneSets.setValue(geneId, genePvalue.get(geneId));

		// Do we need to rank? Rank them by ascending p-value
		GeneSetsRanked geneSetsRanked = null;
		geneSets.setVerbose(verbose);
		if (enrichmentAlgorithmType.isRank()) {
			geneSetsRanked = new GeneSetsRanked(geneSets);
			geneSetsRanked.rankByValue(true);
		}

		//---
		// Run enrichment algorithm
		//---
		EnrichmentAlgorithm algorithm = null;

		switch (enrichmentAlgorithmType) {
		case RANKSUM_GREEDY:
			algorithm = new RankSumPValueGreedyAlgorithm(geneSetsRanked, numberofGeneSetsToSelect);
			break;

		case RANKSUM:
			algorithm = new RankSumPValueAlgorithm(geneSetsRanked, numberofGeneSetsToSelect);
			break;

		case FISHER_GREEDY:
			algorithm = new FisherPValueGreedyAlgorithm(geneSets, numberofGeneSetsToSelect);
			break;

		case FISHER:
			algorithm = new FisherPValueAlgorithm(geneSets, numberofGeneSetsToSelect);
			break;

		case LEADING_EDGE_FRACTION:
			algorithm = new LeadingEdgeFractionAlgorithm(geneSets, numberofGeneSetsToSelect);
			break;

		default:
			throw new RuntimeException("Unimplemented algorithm!");
		}

		// Create 'interesting' genes
		if (enrichmentAlgorithmType.isBinary()) createInterestingGenes();

		// Initialize algorithm parameters
		algorithm.setMaxGeneSetSize(maxGeneSetSize);
		algorithm.setMinGeneSetSize(minGeneSetSize);
		algorithm.setMaxPValue(maxPvalue);
		algorithm.setVerbose(verbose);
		algorithm.setDebug(debug);
		if (enrichmentAlgorithmType.isRank() && enrichmentAlgorithmType.isGreedy()) ((EnrichmentAlgorithmGreedyVariableSize) algorithm).setInitialSize(initGeneSetSize);

		// Run algorithm
		algorithm.select();
	}

	/**
	 * Initialize: read config, database, etc.
	 */
	void initialize() {
		// Read config file
		if (verbose) Timer.showStdErr("Reading configuration file '" + configFile + "'");
		config = new Config(genomeVer, configFile); // Read configuration
		if (verbose) Timer.showStdErr("done");

		// Read database (if gene level p-values are provided, we don't neet to map p_values to genes (we can skip this step)
		if (genePvalueFile.isEmpty()) {
			if (verbose) Timer.showStdErr("Reading database for genome version '" + genomeVer + "' from file '" + config.getFileSnpEffectPredictor() + "' (this might take a while)");
			config.loadSnpEffectPredictor();
			snpEffectPredictor = config.getSnpEffectPredictor();
			genome = config.getGenome();
			if (verbose) Timer.showStdErr("done");

			// Set upstream-downstream interval length
			snpEffectPredictor.setUpDownStreamLength(upDownStreamLength);

			// Build tree
			if (verbose) Timer.showStdErr("Building interval forest");
			snpEffectPredictor.buildForest();
			if (verbose) Timer.showStdErr("done.");

			// Show some genome stats. Chromosome names are shown, a lot of people has problems with the correct chromosome names.
			if (verbose) Timer.showStdErr("Genome stats :\n" + config.getGenome());
		}

		// Read gene set database
		if (verbose) Timer.showStdErr("Reading MSigDb from file: '" + msigdb + "'");
		geneSets = enrichmentAlgorithmType.isRank() ? new GeneSetsRanked(msigdb) : new GeneSets(msigdb);
		if (verbose) Timer.showStdErr("Done. Total:\n\t" + geneSets.getGeneSetCount() + " gene sets\n\t" + geneSets.getGeneCount() + " genes");
	}

	/**
	 * Map <chr,pos,pValue> to gene
	 */
	void mapToGenes() {
		if (verbose) Timer.showStdErr("Mapping p-values to genes.");

		// Create an auto-hash
		genePvalues = new AutoHashMap<String, PvalueList>(new PvalueList());

		//---
		// Map every chr:pos
		//---
		int unmapped = 0, mappedMultiple = 0;
		for (int i = 0; i < chrPosPvalueList.size(); i++) {
			List<String> geneIds = mapToGenes(chrPosPvalueList.getChromosomeName(i), chrPosPvalueList.getStart(i), chrPosPvalueList.getEnd(i));

			// Update counters
			if (geneIds == null || geneIds.isEmpty()) {
				unmapped++;
				continue; // Nothing to do...
			} else if (geneIds.size() > 1) mappedMultiple++;

			// Add pValue to every geneId
			double pvalue = chrPosPvalueList.getPvalue(i);
			for (String geneId : geneIds) {
				PvalueList gpl = genePvalues.getOrCreate(geneId);
				gpl.setGeneId(geneId);
				gpl.add(pvalue);
			}
		}

		//---
		// Show a summary
		//---
		if (verbose) Timer.showStdErr("Done:" //
				+ "\n\tNumber of p-values       : " + chrPosPvalueList.size() //
				+ "\n\tUnmapped                 : " + unmapped //
				+ "\n\tMapped to multiple genes : " + mappedMultiple //
		);

		if (debug) {
			System.err.println("Mapping Gene to pValue:");
			ArrayList<String> geneIds = new ArrayList<String>(genePvalues.keySet());
			Collections.sort(geneIds);
			for (String geneId : geneIds)
				System.err.println("\t" + genePvalues.get(geneId));
		}
	}

	/**
	 * Map a position to a geneId
	 * @param chr
	 * @param start
	 * @return
	 */
	List<String> mapToGenes(String chr, int start, int end) {
		LinkedList<String> geneIds = new LinkedList<String>();

		// Query 
		Marker m = new Marker(genome.getChromosome(chr), start, end, 1, "");

		// Map only to closest gene?
		if (useClosestGene) {
			Gene gene = snpEffectPredictor.queryClosestGene(m);
			if (gene != null) geneIds.add(useGeneId ? gene.getId() : gene.getGeneName());
			return geneIds;
		}

		// Add all genes to list
		Markers hits = snpEffectPredictor.query(m);
		for (Marker mm : hits) {
			if (mm instanceof Gene) {
				Gene gene = (Gene) mm;
				geneIds.add(useGeneId ? gene.getId() : gene.getGeneName());
			}
		}

		return geneIds;
	}

	/**
	 * Parse command line arguments
	 */
	@Override
	public void parseArgs(String[] args) {
		if (args.length == 0) usage(null);

		// Parse comamnd line 
		for (int i = 0; i < args.length; i++) {
			String arg = args[i];

			// Is it an option?
			if (isOpt(arg)) {

				if (arg.equals("-i")) {
					// Input format
					if ((i + 1) < args.length) inputFormat = InputFormat.valueOf(args[++i].toUpperCase());
					else usage("Missing input format in command line option '-i'");
				} else if (arg.equals("-info")) {
					// INFO field name
					if ((i + 1) < args.length) infoName = args[++i];
					else usage("Missing value in command line option '-info'");
				} else if (arg.equals("-ud") || arg.equalsIgnoreCase("-upDownStreamLen")) {
					// Up-downstream length
					if ((i + 1) < args.length) upDownStreamLength = Gpr.parseIntSafe(args[++i]);
					else usage("Missing value in command line option '-ud'");
				} else if (arg.equals("-genePvalue")) {
					// Method for p-value scoring (gene level)
					if ((i + 1) < args.length) {
						String method = args[++i].toUpperCase();
						pvalueSummary = PvalueSummary.valueOf(method);
					} else usage("Missing value in command line option '-genePvalue'");
				} else if (arg.equals("-algo")) {
					// Algorithm to use
					if ((i + 1) < args.length) {
						String algo = args[++i].toUpperCase();
						enrichmentAlgorithmType = EnrichmentAlgorithmType.valueOf(algo);
					} else usage("Missing value in command line option '-algo'");
				} else if (arg.equals("-genePvalueFile")) {
					// Algorithm to use
					if ((i + 1) < args.length) genePvalueFile = args[++i];
					else usage("Missing value in command line option '-genePvalueFile'");
				} else if (arg.equals("-minSetSize")) minGeneSetSize = Gpr.parseIntSafe(args[++i]);
				else if (arg.equals("-maxSetSize")) maxGeneSetSize = Gpr.parseIntSafe(args[++i]);
				else if (arg.equals("-initSetSize")) initGeneSetSize = Gpr.parseIntSafe(args[++i]);
				else if (arg.equals("-rand")) randIterations = Gpr.parseIntSafe(args[++i]);
				else if (arg.equals("-interesting")) interestingPerc = Gpr.parseDoubleSafe(args[++i]);
				else if (arg.equals("-mapClosestGene")) useClosestGene = true;
				else if (arg.equals("-geneId")) useGeneId = true;

			} else if (genomeVer.isEmpty()) genomeVer = arg;
			else if (msigdb.isEmpty()) msigdb = arg;
			else if (inputFile.isEmpty()) inputFile = arg;
		}

		//---
		// Sanity checks
		//---
		if ((inputFormat == InputFormat.VCF) && infoName.isEmpty() && genePvalueFile.isEmpty()) usage("Missing '-info' comamnd line option.");

		// Check input file
		if (inputFile.isEmpty()) inputFile = "-"; // Default is STDIN
		if (!Gpr.canRead(inputFile)) fatalError("Cannot read input file '" + inputFile + "'");

		if (msigdb.isEmpty()) fatalError("Missing Gene-Sets file");
		if (!Gpr.canRead(msigdb)) fatalError("Cannot read Gene-Sets file '" + msigdb + "'");

		if (genomeVer.isEmpty() && genePvalueFile.isEmpty()) usage("Missing genome version.");

		if (maxGeneSetSize <= 0) usage("MaxSetSize must be a positive number.");
		if (minGeneSetSize >= maxGeneSetSize) usage("MaxSetSize (" + maxGeneSetSize + ") must larger than MinSetSize (" + minGeneSetSize + ").");

		if ((interestingPerc < 0) || (interestingPerc > 1)) usage("Interesting percentile must be in the [0 , 1.0] range.");
	}

	/**
	 * Read gene-pValue file
	 * Format: "geneId \t p_value \n"
	 * 
	 * @param genePvalueFile
	 */
	void readGenePvalues(String genePvalueFile) {
		if (verbose) Timer.showStdErr("Reading gene p-values file '" + genePvalueFile + "'");

		genePvalue = new HashMap<String, Double>();

		// Read the whole file
		String lines[] = Gpr.readFile(genePvalueFile).split("\n");

		// Parse each line
		double minp = 1.0;
		for (String line : lines) {
			String rec[] = line.split("\\s");
			String geneId = rec[0].trim();
			double pValue = Gpr.parseDoubleSafe(rec[1]);

			if ((pValue > 0) && (pValue <= 1.0)) { // Assume that a p-value of zero is a parsing error
				genePvalue.put(geneId, pValue);
				minp = Math.min(minp, pValue);
			} else if (verbose) System.err.println("\tWarning: Ignoring entry (zero p-value):\t'" + line + "'");
		}

		if (verbose) Timer.showStdErr("Done."//
				+ "\n\tAdded       : " + genePvalue.size() //
				+ "\n\tMin p-value : " + minp //
		);
	}

	/**
	 * Read input file and populate 'chrPosPvalueList'
	 */
	void readInput() {
		if (verbose) Timer.showStdErr("Reading input file '" + inputFile + "' (format '" + inputFormat + "')");

		switch (inputFormat) {
		case VCF:
			chrPosPvalueList = readInputVcf();
			break;

		case TXT:
			chrPosPvalueList = readInputTxt();
			break;

		case BED:
			chrPosPvalueList = readInputBed();
			break;

		default:
			fatalError("Input format '" + inputFormat + "' not supported!");
		}

		if (verbose) {
			System.err.println("");
			Timer.showStdErr("Done.");
		}

		if (debug) {
			// Show data
			System.err.println("P-values:\n\tchr\tstart\tend\tp_value");
			for (int i = 0; i < chrPosPvalueList.size(); i++)
				System.err.println("\t" + chrPosPvalueList.getChromosomeName(i) + "\t" + chrPosPvalueList.getStart(i) + "\t" + chrPosPvalueList.getEnd(i) + "\t" + chrPosPvalueList.getPvalue(i));
		}
	}

	/**
	 * Read input in BED format
	 * 
	 * Format: "chr \t start \t end \t id \t pValue \n"
	 *         start : zero-based
	 *         end   : zero-based open
	 *
	 */
	ChrPosPvalueList readInputBed() {
		ChrPosPvalueList cppList = new ChrPosPvalueList();

		int num = 1;
		BedFileIterator bfi = new BedFileIterator(inputFile);
		for (SeqChange sc : bfi) {
			cppList.add(sc.getChromosome(), sc.getStart(), sc.getEnd(), sc.getScore());
			if (verbose) Gpr.showMark(num++, READ_INPUT_SHOW_EVERY);
		}

		return cppList;
	}

	/**
	 * Read input in TXT format
	 * 
	 * Format: "chr \t pos \t p-value \n"
	 * 
	 * Note: Position is 1-based coordinate
	 */
	ChrPosPvalueList readInputTxt() {
		ChrPosPvalueList cppList = new ChrPosPvalueList();
		Genome genome = config.getGenome();

		int num = 1;
		LineFileIterator lfi = new LineFileIterator(inputFile);
		for (String line : lfi) {
			if (line.startsWith("#")) continue; // Ignore lines that start with '#'
			String fields[] = line.split("\t");

			// Sanity check
			if (fields.length < 4) {
				System.err.println("Warning: Ignoring line number " + lfi.getLineNum() + "." //
						+ " Exepcting format 'chr\tpos\tp_value\n'.\n" //
						+ "\tLine:\t'" + line + "'" //
				);
				continue;
			}

			// Parse fields
			String chr = fields[0];
			int start = Gpr.parseIntSafe(fields[1]) - 1; // Input format is 1-based
			int end = Gpr.parseIntSafe(fields[2]) - 1; // Input format is 1-based
			double pvalue = Gpr.parseDoubleSafe(fields[3]);

			// Add data to list
			Chromosome chromo = genome.getOrCreateChromosome(chr);
			cppList.add(chromo, start, end, pvalue);

			if (verbose) Gpr.showMark(num++, READ_INPUT_SHOW_EVERY);
		}

		return cppList;

	}

	/**
	 * Read input in VCF format
	 */
	ChrPosPvalueList readInputVcf() {
		ChrPosPvalueList cppList = new ChrPosPvalueList();

		int num = 1;
		VcfFileIterator vcf = new VcfFileIterator(inputFile);
		for (VcfEntry ve : vcf) {
			double pvalue = ve.getInfoFloat(infoName);

			if (Double.isNaN(pvalue)) {
				// Error
				System.err.println("Warning: Cannot find INFO field '" + infoName + "'.\n\tIgnoring VCF entry." + vcf.getLineNum() + "\n\t" + ve);
			} else {
				// Add to list
				cppList.add(ve.getChromosome(), ve.getStart(), ve.getEnd(), pvalue);
			}

			if (verbose) Gpr.showMark(num++, READ_INPUT_SHOW_EVERY);
		}

		return cppList;
	}

	/**
	* Run command
	*/
	@Override
	public boolean run() {
		initialize();

		if (genePvalueFile.isEmpty()) {
			// Perform 'normal' procedure
			readInput(); // Read input files (p-values)
			mapToGenes(); // Map <chr,pos,pValue> to gene
			scoreGenes(); // Get one score (pValue) per gene
			correctPvalues(); // Correct gene scores
		} else {
			// P-values already mapped to genes, provided in a file
			readGenePvalues(genePvalueFile);
		}

		enrichmentAnalysis(); // Perform enrichment analysis
		if (randIterations > 0) runRand(); // Perform random iterations

		if (verbose) Timer.showStdErr("Done.");
		return true;
	}

	/**
	* Run enrichment analysis using random p-values
	*/
	public boolean runRand() {
		HashMap<String, Double> genePvalueOri = genePvalue; // Save original p-values

		for (int iter = 1; iter <= randIterations; iter++) {
			Timer.showStdErr("Random p-values. Iteration " + iter);

			// Create random pvalues based on input 
			genePvalue = new HashMap<String, Double>();
			for (String gene : genePvalueOri.keySet())
				genePvalue.put(gene, Math.random());

			// Perform enrichment analysis
			enrichmentAnalysis();
		}

		genePvalue = genePvalueOri; // Restore original values
		if (verbose) Timer.showStdErr("Done.");
		return true;
	}

	/**
	 * Get one score (pValue) per gene
	 */
	void scoreGenes() {
		if (verbose) Timer.showStdErr("Aggregating p-values by gene (scoring genes)");

		// Create one pValue per gene
		genePvalue = new HashMap<String, Double>();
		if (debug) System.err.println("\tp-value\tgeneId");
		for (String geneId : genePvalues.keySet()) {
			// Calculate aggregated score
			PvalueList gpl = genePvalues.get(geneId);
			double pValue = gpl.pValue(pvalueSummary);

			// Add to map
			genePvalue.put(geneId, pValue);
			if (debug) System.err.println(String.format("\t%.2e\t%s", pValue, geneId));
		}

		if (verbose) Timer.showStdErr("Done.");
	}

	/**
	 * Show usage and exit
	 */
	@Override
	public void usage(String message) {
		if (message != null) System.err.println("Error: " + message + "\n");
		System.err.println("snpEff version " + SnpEff.VERSION);
		System.err.println("Usage: snpEff gsa [options] genome_version geneSets.gmt input_file");
		System.err.println("\n\tInput data options:");
		System.err.println("\t-geneId            : Use geneID instead of gene names. Default: " + useGeneId);
		System.err.println("\t-i <format>        : Input format {vcf, bed, txt}. Default: " + inputFormat);
		System.err.println("\t-info <name>       : INFO tag used for p-values (in VCF input format).");
		System.err.println("\n\tAlgorithm options:");
		System.err.println("\t-algo <name>       : Gene set enrichment algorithm {FISHER_GREEDY, RANKSUM_GREEDY, FISHER, RANKSUM}. Default: " + enrichmentAlgorithmType);
		System.err.println("\t-genePvalue        : Method to summarize gene p-values {MIN, AVG, AVG10, FISHER_CHI_SQUARE, Z_SCORES, SIMES}. Default: " + pvalueSummary);
		System.err.println("\t-genePvalueCorr    : Correction method for gene-summarized p-values {NONE}. Default: " + correctionMethod);
		System.err.println("\t-mapClosestGene    : Map to closest gene. Default: " + useClosestGene);
		System.err.println("\t-rand <num>        : Perform 'num' iterations using random p-values. Default: " + randIterations);
		System.err.println("\n\tAlgorithm specific options: FISHER and FISHER_GREEDY");
		System.err.println("\t-interesting <num> : Consider a gene 'interesting' if the p-value is in the 'num' percentile. Default: " + interestingPerc);
		System.err.println("\n\tGene Set options:");
		System.err.println("\t-minSetSize <num>  : Minimum number of genes in a gene set. Default: " + minGeneSetSize);
		System.err.println("\t-maxSetSize <num>  : Maximum number of genes in a gene set. Default: " + maxGeneSetSize);
		System.err.println("\t-initSetSize <num> : Initial number of genes in a gene set (size range algorithm). Default: " + initGeneSetSize);
		System.exit(-1);
	}
}
