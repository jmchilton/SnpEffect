package ca.mcgill.mcb.pcingola.snpEffect.commandLine;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import net.sf.samtools.AbstractBAMFileIndex;
import net.sf.samtools.BAMIndexMetaData;
import net.sf.samtools.SAMFileReader;
import ca.mcgill.mcb.pcingola.coverage.CountReadsOnMarkers;
import ca.mcgill.mcb.pcingola.fileIterator.BedFileIterator;
import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.SnpEffectPredictor;
import ca.mcgill.mcb.pcingola.stats.ReadsOnMarkersModel;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * Count reads from a BAM file given a list of intervals
 * 
 * @author pablocingolani
 */
public class SnpEffCmdCountReads extends SnpEff {

	public static boolean debug = true;

	CountReadsOnMarkers countReadsOnMarkers;
	ReadsOnMarkersModel readsOnMarkersModel;
	SnpEffectPredictor snpEffectPredictor;
	List<String> samFileNames; // SAM or BAM files 
	List<String> bedFileNames; // BED files for custom intervals

	public SnpEffCmdCountReads() {
		samFileNames = new ArrayList<String>();
		bedFileNames = new ArrayList<String>();
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
			if (args[i].equals("-i")) bedFileNames.add(args[++i]);
			else if ((genomeVer == null) || genomeVer.isEmpty()) genomeVer = args[i];
			else samFileNames.add(args[i]);
		}

		// Sanity check
		if ((genomeVer == null) || genomeVer.isEmpty()) usage("Missing genome version");
		if (samFileNames.size() < 1) usage("Missing SAM/BAM file/s");
	}

	/**
	 * Calculate pvalues for 
	 */
	void pvalues() {
		readsOnMarkersModel = new ReadsOnMarkersModel(snpEffectPredictor);
		int readLength = countReadsOnMarkers.getReadLengthAvg();
		if (verbose) Timer.showStdErr("Calculating probability model for read length " + readLength);
		readsOnMarkersModel.setReadLength(readLength);
		readsOnMarkersModel.setVerbose(verbose);
		readsOnMarkersModel.setMarkerTypes(countReadsOnMarkers.getMarkerTypes());

		// Run model
		readsOnMarkersModel.run();
		System.err.println(countReadsOnMarkers.probabilityTable(readsOnMarkersModel.getProb()));
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
		snpEffectPredictor = config.getSnpEffectPredictor(); // Read snpEffect predictor
		countReadsOnMarkers = new CountReadsOnMarkers(snpEffectPredictor);
		if (verbose) Timer.showStdErr("done");

		// Load BED files
		for (String bedFile : bedFileNames) {
			// Load file
			if (verbose) Timer.showStdErr("Reading intervals from file '" + bedFile + "'");
			String baseName = Gpr.removeExt(Gpr.baseName(bedFile));

			BedFileIterator bed = new BedFileIterator(bedFile);
			List<SeqChange> markers = bed.load();
			for (SeqChange marker : markers) {
				marker.setId(baseName + ":" + marker.getId());
				snpEffectPredictor.add(marker);
				countReadsOnMarkers.addMarkerType(marker, baseName);
			}
			if (verbose) Timer.showStdErr("Done. Intervals added : " + markers.size());
		}

		// Build forest
		if (verbose) Timer.showStdErr("Building interval forest");
		snpEffectPredictor.buildForest();
		if (verbose) Timer.showStdErr("done");

		// Count reads
		countReadsOnMarkers.setVerbose(verbose);
		for (String file : samFileNames)
			countReadsOnMarkers.addFile(file);
		countReadsOnMarkers.count();

		if (!quiet) {
			countReadsOnMarkers.print();
			pvalues(); // Calculate (and show) pvalues
		}

		return true;
	}

	@Override
	public void usage(String message) {
		if (message != null) System.err.println("Error: " + message + "\n");
		System.err.println("snpEff version " + VERSION);
		System.err.println("Usage: snpEff countReads [-i intervals.bed] genome readsFile_1 readsFile_2 ...  readsFile_N");
		System.err.println("\t-i intervals.bed : User defined intervals. Mutiple '-i' commands are allowed.");
		System.err.println("\treadsFile        : A file contianing the reads. Either BAM or SAM format.");
		System.exit(-1);
	}
}
