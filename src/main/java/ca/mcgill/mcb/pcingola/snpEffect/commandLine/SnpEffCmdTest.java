package ca.mcgill.mcb.pcingola.snpEffect.commandLine;

import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Markers;
import ca.mcgill.mcb.pcingola.interval.Motif;
import ca.mcgill.mcb.pcingola.probablility.FisherExactTest;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.SnpEffectPredictor;
import ca.mcgill.mcb.pcingola.stats.CountByType;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * Command line: Test
 * 
 * Note: Used for testing weird stuff
 * 
 * @author pcingola
 */
public class SnpEffCmdTest extends SnpEff {

	public static final double P_VALUE_LIMIT = 0.05;

	SnpEffectPredictorLoader sepLoader;
	SnpEffectPredictor snpEffectPredictor;
	String inputFile;

	public SnpEffCmdTest() {
		super();
		sepLoader = new SnpEffectPredictorLoader();
	}

	/**
	 * Parse command line arguments
	 */
	@Override
	public void parseArgs(String[] args) {
		if (args.length < 2) usage(null);

		// Parse command line args in SnpEffectPredictorLoader
		sepLoader.setVerbose(verbose);
		sepLoader.setDebug(debug);
		sepLoader.setQuiet(quiet);
		args = sepLoader.parseArgs(args);

		// Parse all other command line options
		for (int idx = 0; idx < args.length; idx++) {
			String arg = args[idx];

			if (isOpt(arg)) {
				usage("Unknown opton '" + arg + "'");
			} else if (inputFile == null) inputFile = args[idx];
		}

		// Sanity check
		if (inputFile == null) usage("Missing input file");
	}

	/**
	 * Print a counter using a label on each line
	 * @param label
	 * @param countByType
	 */
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
		// Initialize
		//---
		if (verbose) Timer.showStdErr("Reading configuration file '" + configFile + "'");
		config = new Config(sepLoader.getGenomeVer(), configFile); // Read configuration
		if (verbose) Timer.showStdErr("done");

		sepLoader.load(config); // Load database, build forest

		//---
		// Load intervals
		//---
		if (verbose) Timer.showStdErr("Loading intervals from '" + inputFile + "'");
		Markers intervals = Markers.readMarkers(inputFile);
		if (verbose) Timer.showStdErr("Done. Intervals added : " + intervals.size());

		//---
		// Process input file
		//---
		sepLoader.build();
		run(intervals);

		if (verbose) Timer.showStdErr("Done");
		return true;
	}

	void run(Markers intervals) {
		if (verbose) Timer.showStdErr("Runnig");
		SnpEffectPredictor sep = config.getSnpEffectPredictor();

		//---
		// Count for each motif
		//---
		CountByType countByMotif = new CountByType();

		for (Marker m : intervals) {
			System.out.println(m);

			Markers results = sep.query(m);
			for (Marker r : results) {
				if (r instanceof Motif) {
					Motif motif = (Motif) r;
					System.out.println("\t" + motif.toString() + "\t" + motif.getPwmId() + "\t" + motif.getPwmName());
					countByMotif.inc(motif.getPwmId());
				}
			}
		}

		//---
		// Show results
		//---

		// Count totals by motif
		CountByType totalByMotif = new CountByType();
		for (Marker m : sep.getMarkers())
			if (m instanceof Motif) {
				Motif motif = (Motif) m;
				totalByMotif.inc(motif.getPwmId());
			}

		// Count motifs in intervals
		long totalHits = countByMotif.sum();
		long total = totalByMotif.sum();
		System.out.println("Total hits        : " + totalHits);
		System.out.println("Total motif sites : " + total);

		// Calculate p-values and show each motif
		for (String id : countByMotif.keysSorted()) {
			long count = countByMotif.get(id);

			int k = (int) count;
			int N = (int) total;
			int D = (int) totalByMotif.get(id);
			int n = (int) totalHits;
			double pvalue = FisherExactTest.get().pValueUp(k, N, D, n, P_VALUE_LIMIT);

			System.out.println(id + "\t" + count + "\t" + totalByMotif.get(id) + "\t" + pvalue);
		}
	}

	/**
	 * Show usage and exit
	 */
	@Override
	public void usage(String message) {
		if (message != null) System.err.println("Error: " + message + "\n");
		System.err.println("snpEff version " + SnpEff.VERSION);
		System.err.println("Usage: snpEff test [options] genomeVer inputFile");
		System.err.println("Where : ");
		System.err.println("\tinputFile : Input intervals file. Can be BED, TXT, BigBed or VCF");
		System.err.println("Options : ");
		sepLoader.usage(null); // Show SnpEffect predictor loader options

		System.exit(-1);
	}

}
