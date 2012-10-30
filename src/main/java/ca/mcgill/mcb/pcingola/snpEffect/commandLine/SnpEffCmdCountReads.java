package ca.mcgill.mcb.pcingola.snpEffect.commandLine;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import net.sf.samtools.AbstractBAMFileIndex;
import net.sf.samtools.BAMIndexMetaData;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Genome;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.SnpEffectPredictor;
import ca.mcgill.mcb.pcingola.stats.CountByType;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * Count reads from a BAM file given a list of intervals
 * 
 * @author pablocingolani
 */
public class SnpEffCmdCountReads extends SnpEff {

	public static int SHOW_EVERY = 10;
	public static boolean debug = true;

	List<String> samFileNames;
	ArrayList<CountByType> countByFile;
	SnpEffectPredictor snpEffectPredictor;
	boolean verbose = false;

	public SnpEffCmdCountReads() {
		samFileNames = new ArrayList<String>();
		countByFile = new ArrayList<CountByType>();
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
	 * A list of all IDs and parent IDs until chromosome
	 * @param m
	 * @return
	 */
	String idChain(Marker m) {
		StringBuilder sb = new StringBuilder();
		for (; (m != null) && !(m instanceof Chromosome); m = m.getParent()) {
			if (sb.length() > 0) sb.append(";");

			switch (m.getType()) {
			case EXON:
				sb.append("exon_" + ((Exon) m).getRank());
				break;

			case TRANSCRIPT:
			case GENE:
				sb.append(m.getId());
				break;

			default:
				sb.append(m.getClass().getSimpleName());
				break;
			}
		}
		return sb.toString();
	}

	/**
	 * Parse
	 * @param args
	 */
	@Override
	public void parseArgs(String[] args) {
		// Parse command line arguments
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-v")) verbose = true;
			else if ((genomeVer == null) || genomeVer.isEmpty()) genomeVer = args[i];
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

		runCountIntervals();
		System.out.print(this);

		return true;
	}

	/**
	 * Count reads onto intervals
	 */
	void runCountIntervals() {
		Genome genome = config.getGenome();

		// Iterate over all BAM/SAM files
		for (String samFileName : samFileNames) {
			try {
				if (verbose) Timer.showStdErr("Reading reads file '" + samFileName + "'");
				CountByType countReads = new CountByType();

				// Open file
				SAMFileReader sam = new SAMFileReader(new File(samFileName));
				for (SAMRecord samRecord : sam) {
					System.out.println(samRecord + "\t\t" + samRecord.getReferenceName() + ":" + samRecord.getAlignmentStart() + "-" + samRecord.getAlignmentEnd());

					if (!samRecord.getReadUnmappedFlag()) { // Mapped?
						Chromosome chr = genome.getChromosome(samRecord.getReferenceName());
						if (chr != null) {
							// Create a marker
							Marker read = new Marker(chr, samRecord.getAlignmentStart(), samRecord.getAlignmentEnd(), 1, "");

							// Find all intersects
							Set<Marker> regions = snpEffectPredictor.regionsMarkers(read);
							for (Marker s : regions) {
								String idChain = idChain(s);
								countReads.inc(idChain);
								System.out.println("\t\t" + idChain);
							}
						}
					}
				}
				sam.close();

				countByFile.add(countReads); // Add count to list
			} catch (Exception e) {
				e.printStackTrace();
			}

			if (verbose) {
				System.err.println("");
				Timer.showStdErr("Finished reding file " + samFileName);
			}
		}
		if (verbose) Timer.showStdErr("Done.");
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();

		// Show title
		sb.append("chr\tstart\tend");
		for (int j = 0; j < countByFile.size(); j++)
			sb.append("\t" + samFileNames.get(j));
		sb.append("\n");

		return sb.toString();
	}

	@Override
	public void usage(String message) {
		if (message != null) System.err.println("Error: " + message + "\n");
		System.err.println("snpEff version " + VERSION);
		System.err.println("Usage: snpEff countReads [-v] genome readsFile_1 readsFile_2 ...  readsFile_N");
		System.err.println("\treadsFile : A file contianing the reads. Either BAM or SAM format.");
		System.exit(-1);
	}
}
