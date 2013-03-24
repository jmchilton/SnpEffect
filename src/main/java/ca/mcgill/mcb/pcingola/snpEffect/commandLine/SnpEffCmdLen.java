package ca.mcgill.mcb.pcingola.snpEffect.commandLine;

import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.SnpEffectPredictor;
import ca.mcgill.mcb.pcingola.stats.ReadsOnMarkersModel;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * Calculate the maximum interval length by type, for all markers in a genome
 * 
 * 
 * @author pcingola
 */
public class SnpEffCmdLen extends SnpEff {

	int readLength, numIterations, numReads;
	ReadsOnMarkersModel readsOnMarkersModel;

	public SnpEffCmdLen() {
		super();
	}

	@Override
	public void parseArgs(String[] args) {
		this.args = args;
		for (int i = 0; i < args.length; i++) {

			// Argument starts with '-'?
			if (args[i].equals("-r")) {
				if ((i + 1) < args.length) readLength = Gpr.parseIntSafe(args[++i]);
				else usage("Missing value for parameter '-r'");
			} else if (args[i].equals("-iter")) {
				if ((i + 1) < args.length) numIterations = Gpr.parseIntSafe(args[++i]);
				else usage("Missing value for parameter '-iter'");
			} else if (args[i].equals("-reads")) {
				if ((i + 1) < args.length) numReads = Gpr.parseIntSafe(args[++i]);
				else usage("Missing value for parameter '-reads'");
			} else if (args[i].equals("-r")) {
				if ((i + 1) < args.length) readLength = Gpr.parseIntSafe(args[++i]);
				else usage("Missing value for parameter '-r'");

			} else if (genomeVer.isEmpty()) genomeVer = args[i];
			else usage("Unknow parameter '" + args[i] + "'");
		}

		// Check: Do we have all required parameters?
		if (genomeVer.isEmpty()) usage("Missing genomer_version parameter");
		if (readLength < 0) usage("Read length should be a non-negative number");
	}

	/**
	 * Run
	 * @return 
	 */
	@Override
	public boolean run() {
		if (verbose) Timer.showStdErr("Loading config");
		Config config = new Config(genomeVer);

		if (verbose) Timer.showStdErr("Loading predictor");
		SnpEffectPredictor snpEffectPredictor = config.loadSnpEffectPredictor();

		if (verbose) Timer.showStdErr("Building interval forest");
		snpEffectPredictor.buildForest();

		readsOnMarkersModel = new ReadsOnMarkersModel(snpEffectPredictor);
		readsOnMarkersModel.setVerbose(verbose);

		if (verbose) Timer.showStdErr("Counting bases");
		readsOnMarkersModel.run(); // Count 
		if (!quiet) System.out.println(readsOnMarkersModel);

		// Perform some random sampling
		if ((numIterations > 0) && (readLength > 0)) readsOnMarkersModel.randomSampling(numIterations, readLength, numReads);

		return true;
	}

	/**
	 * Show usage and exit
	 */
	@Override
	public void usage(String message) {
		if (message != null) System.err.println("Error: " + message + "\n");
		System.err.println("snpEff version " + SnpEff.VERSION);
		System.err.println("Usage: snpEff len [options] genome_version");
		System.err.println("Options:");
		System.err.println("\t-r     <num> : Assume a read size of 'num' bases.");
		System.err.println("\t-iter  <num> : Perform 'num' iterations of random sampling.");
		System.err.println("\t-reads <num> : Each random sampling iteration has 'num' reads.");
		System.exit(-1);
	}

}
