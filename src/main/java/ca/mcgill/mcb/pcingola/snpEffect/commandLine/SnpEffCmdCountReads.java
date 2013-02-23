package ca.mcgill.mcb.pcingola.snpEffect.commandLine;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import net.sf.samtools.AbstractBAMFileIndex;
import net.sf.samtools.BAMIndexMetaData;
import net.sf.samtools.SAMFileReader;
import ca.mcgill.mcb.pcingola.coverage.CountReadsOnMarkers;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.SnpEffectPredictor;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * Count reads from a BAM file given a list of intervals
 * 
 * @author pablocingolani
 */
public class SnpEffCmdCountReads extends SnpEff {

	public static boolean debug = true;

	CountReadsOnMarkers countReadsOnMarkers;
	SnpEffectPredictor snpEffectPredictor;
	List<String> samFileNames;

	public SnpEffCmdCountReads() {
		samFileNames = new ArrayList<String>();
	}

	/**
	 * Count all reads in a BAM file 
	 * Note: It uses the BAM index
	 * 
	 * @param samReader
	 * @return
	 */
	int countTotalReads(String samFileName) {
		try {
			if (verbose) Timer.showStdErr("Counting reads on file: " + samFileName);
			SAMFileReader samReader = new SAMFileReader(new File(samFileName));
			AbstractBAMFileIndex index = (AbstractBAMFileIndex) samReader.getIndex();
			int count = 0;
			for (int i = 0; i < index.getNumberOfReferences(); i++) {
				BAMIndexMetaData meta = index.getMetaData(i);
				count += meta.getAlignedRecordCount();
			}
			samReader.close();
			if (verbose) Timer.showStdErr("Total " + count + " reads.");
			return count;
		} catch (Exception e) {
			// Error? (e.g. no index)
			System.err.println("ERROR! BAM file not indexed?");
			return -1;
		}
	}

	/**
	 * Parse
	 * @param args
	 */
	@Override
	public void parseArgs(String[] args) {
		// Parse command line arguments
		for (int i = 0; i < args.length; i++) {
			if ((genomeVer == null) || genomeVer.isEmpty()) genomeVer = args[i];
			else samFileNames.add(args[i]);
		}

		// Sanity check
		if ((genomeVer == null) || genomeVer.isEmpty()) usage("Missing genome version");
		if (samFileNames.size() < 1) usage("Missing SAM/BAM file/s");
	}

	/**
	 * Run
	 */
	@Override
	public boolean run() {
		// Load Config
		if (verbose) Timer.showStdErr("Reading configuration file '" + configFile + "'");
		config = new Config(genomeVer, configFile); // Read configuration
		if (verbose) Timer.showStdErr("done");

		// Load database
		if (verbose) Timer.showStdErr("Reading database for genome '" + genomeVer + "'");
		config.loadSnpEffectPredictor(); // Read snpEffect predictor
		if (verbose) Timer.showStdErr("done");

		// Build forest
		if (verbose) Timer.showStdErr("Building interval forest");
		snpEffectPredictor = config.getSnpEffectPredictor(); // Read snpEffect predictor
		snpEffectPredictor.buildForest();
		if (verbose) Timer.showStdErr("done");

		// Count reads
		CountReadsOnMarkers countReadsOnMarkers = new CountReadsOnMarkers(snpEffectPredictor);
		countReadsOnMarkers.setVerbose(verbose);
		for (String file : samFileNames)
			countReadsOnMarkers.addFile(file);
		countReadsOnMarkers.count();
		countReadsOnMarkers.print();

		return true;
	}

	@Override
	public void usage(String message) {
		if (message != null) System.err.println("Error: " + message + "\n");
		System.err.println("snpEff version " + VERSION);
		System.err.println("Usage: snpEff countReads genome readsFile_1 readsFile_2 ...  readsFile_N");
		System.err.println("\treadsFile : A file contianing the reads. Either BAM or SAM format.");
		System.exit(-1);
	}
}
