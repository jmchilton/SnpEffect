package ca.mcgill.mcb.pcingola.coverage;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import ca.mcgill.mcb.pcingola.fileIterator.BedFileIterator;
import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.interval.Genome;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Markers;
import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.outputFormatter.OutputFormatter;
import ca.mcgill.mcb.pcingola.probablility.Binomial;
import ca.mcgill.mcb.pcingola.snpEffect.SnpEffectPredictor;
import ca.mcgill.mcb.pcingola.stats.CountByKey;
import ca.mcgill.mcb.pcingola.stats.CountByType;
import ca.mcgill.mcb.pcingola.stats.CoverageByType;
import ca.mcgill.mcb.pcingola.stats.PosStats;
import ca.mcgill.mcb.pcingola.stats.plot.GoogleGeneRegionChart;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;

/**
 * Count how many reads map (from many SAM/BAM files) onto markers
 * @author pcingola
 */
public class CountReadsOnMarkers {

	public static int SHOW_EVERY = 10000;
	public static boolean debug = true;

	boolean verbose = false; // Be verbose
	int countExceptions = 0;
	int readLengthCount;
	long readLengthSum;
	List<String> samFileNames;
	List<String> names;
	ArrayList<CountByKey<Marker>> countReadsByFile;
	ArrayList<CountByKey<Marker>> countBasesByFile;
	ArrayList<CountByType> countTypesByFile;
	SnpEffectPredictor snpEffectPredictor;
	HashMap<Marker, String> markerTypes;
	ArrayList<CoverageByType> coverageByFile;
	Genome genome;

	public CountReadsOnMarkers() {
		init(null);
	}

	public CountReadsOnMarkers(SnpEffectPredictor snpEffectPredictor) {
		init(snpEffectPredictor);
	}

	/**
	 * Add a SAM/BAM file to be processed
	 * @param samFileName
	 */
	public void addFile(String samFileName) {
		samFileNames.add(samFileName);
		names.add(Gpr.removeExt(Gpr.baseName(samFileName)));
	}

	public void addMarkerType(Marker marker, String type) {
		markerTypes.put(marker, type);
	}

	/**
	 * Create a collection of all markers
	 * @return
	 */
	List<Marker> allMarkers() {
		// Retrieve all possible keys, sort them
		HashSet<Marker> keys = new HashSet<Marker>();
		for (CountByKey<Marker> cbt : countReadsByFile)
			keys.addAll(cbt.keySet());

		ArrayList<Marker> keysSorted = new ArrayList<Marker>(keys.size());
		keysSorted.addAll(keys);
		Collections.sort(keysSorted);
		return keysSorted;
	}

	/**
	 * Count markers from all files
	 */
	public void count() {
		genome = snpEffectPredictor.getGenome();

		readLengthSum = 0;
		readLengthCount = 0;

		// Iterate over all BAM/SAM files
		for (String samFileName : samFileNames) {
			try {
				if (verbose) Timer.showStdErr("Reading reads file '" + samFileName + "'");
				CountByKey<Marker> countReads = new CountByKey<Marker>();
				CountByKey<Marker> countBases = new CountByKey<Marker>();
				CountByType countTypes = new CountByType();
				CoverageByType coverageByType = new CoverageByType();

				countFile(samFileName, countReads, countBases, countTypes, coverageByType);

				// Add count to list
				countReadsByFile.add(countReads);
				countBasesByFile.add(countBases);
				countTypesByFile.add(countTypes);
				coverageByFile.add(coverageByType);
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

	/**
	 * Count all markers from a BED file
	 */
	void countBedFile(String fileName, CountByKey<Marker> countReads, CountByKey<Marker> countBases, CountByType countTypes, CoverageByType coverageByType) {
		// Open file
		int readNum = 1;

		for (SeqChange read : new BedFileIterator(fileName)) {
			try {
				readLengthCount++;
				readLengthSum += read.size();
				countMarker(read, countReads, countBases, countTypes, coverageByType);
				if (verbose) Gpr.showMark(readNum, SHOW_EVERY);
				readNum++;
			} catch (Exception e) {
				countExceptions++;
				if (countExceptions < 10) e.printStackTrace();
				else if (countExceptions == 10) System.err.println("Not showing more exceptions!");
			}
		}
	}

	/**
	 * Count all markers from a SAM/BAM file
	 */
	void countFile(String fileName, CountByKey<Marker> countReads, CountByKey<Marker> countBases, CountByType countTypes, CoverageByType coverageByType) {
		String fl = fileName.toLowerCase();

		if (fl.endsWith(".bam") || fl.endsWith(".sam")) countSamFile(fileName, countReads, countBases, countTypes, coverageByType);
		else if (fl.endsWith(".vcf") || fl.endsWith(".vcf.gz")) countVcfFile(fileName, countReads, countBases, countTypes, coverageByType);
		else if (fl.endsWith(".bed") || fl.endsWith(".bed.gz")) countBedFile(fileName, countReads, countBases, countTypes, coverageByType);
		else throw new RuntimeException("Unrecognized file extention. Supported types: BAM, SAM, BED, VCF.");
	}

	/**
	 * Count one marker
	 */
	void countMarker(Marker read, CountByKey<Marker> countReads, CountByKey<Marker> countBases, CountByType countTypes, CoverageByType coverageByType) {
		// Find all intersects
		Markers regions = snpEffectPredictor.queryDeep(read);

		HashSet<String> doneClass = new HashSet<String>();
		for (Marker m : regions) {
			countReads.inc(m); // Count reads
			countBases.inc(m, m.intersectSize(read)); // Count number bases that intersect

			// Count by marker type (make sure we only count once per read)
			String type = markerType(m);
			if (!doneClass.contains(type)) {
				countTypes.inc(type); // Count reads
				doneClass.add(type); // Do not count twice

				PosStats posStats = coverageByType.getOrCreate(type);
				posStats.sample(read, m);
			}

		}

	}

	/**
	 * Count how many of each marker type are there
	 * @return
	 */
	CountByType countMarkerTypes(Collection<Marker> markersToCount) {
		CountByType countByMarkerType = new CountByType();
		for (Marker marker : markersToCount)
			countByMarkerType.inc(markerType(marker));
		return countByMarkerType;
	}

	/**
	 * Count all markers from a SAM/BAM file
	 */
	void countSamFile(String fileName, CountByKey<Marker> countReads, CountByKey<Marker> countBases, CountByType countTypes, CoverageByType coverageByType) {
		// Open file
		int readNum = 1;
		SAMFileReader sam = new SAMFileReader(new File(fileName));
		sam.setValidationStringency(ValidationStringency.SILENT);

		for (SAMRecord samRecord : sam) {
			try {
				if (!samRecord.getReadUnmappedFlag()) { // Mapped?
					Chromosome chr = genome.getOrCreateChromosome(samRecord.getReferenceName());
					if (chr != null) {
						// Create a marker from read
						Marker read = new Marker(chr, samRecord.getAlignmentStart(), samRecord.getAlignmentEnd(), 1, "");
						readLengthCount++;
						readLengthSum += read.size();

						countMarker(read, countReads, countBases, countTypes, coverageByType);
					}
				}

				if (verbose) Gpr.showMark(readNum, SHOW_EVERY);
				readNum++;
			} catch (Exception e) {
				countExceptions++;
				if (countExceptions < 10) e.printStackTrace();
				else if (countExceptions == 10) System.err.println("Not showing more exceptions!");
			}

		}
		sam.close();
	}

	/**
	 * Count all markers from a VCF file
	 */
	void countVcfFile(String fileName, CountByKey<Marker> countReads, CountByKey<Marker> countBases, CountByType countTypes, CoverageByType coverageByType) {
		// Open file
		int readNum = 1;

		for (VcfEntry read : new VcfFileIterator(fileName)) {
			try {
				readLengthCount++;
				readLengthSum += read.size();
				countMarker(read, countReads, countBases, countTypes, coverageByType);
				if (verbose) Gpr.showMark(readNum, SHOW_EVERY);
				readNum++;
			} catch (Exception e) {
				countExceptions++;
				if (countExceptions < 10) e.printStackTrace();
				else if (countExceptions == 10) System.err.println("Not showing more exceptions!");
			}
		}
	}

	public HashMap<Marker, String> getMarkerTypes() {
		return markerTypes;
	}

	/**
	 * Average read length
	 * @return
	 */
	public int getReadLengthAvg() {
		if (readLengthCount <= 0) return 0;
		double rl = ((double) readLengthSum) / readLengthCount;
		return (int) Math.round(rl);
	}

	public String html() {
		StringBuilder sbHead = new StringBuilder();
		StringBuilder sbBody = new StringBuilder();

		for (int i = 0; i < names.size(); i++) {
			String name = names.get(i);
			Gpr.debug("File : " + i + "\t" + name);
			CoverageByType cvt = coverageByFile.get(i);

			GoogleGeneRegionChart grc = new GoogleGeneRegionChart(cvt, name);
			sbHead.append(grc.toStringHtmlHeader());
			sbBody.append(grc.toStringHtmlBody());
		}

		return sbHead.toString() + sbBody.toString();
	}

	/**
	 * Initialize
	 * @param snpEffectPredictor
	 */
	void init(SnpEffectPredictor snpEffectPredictor) {
		samFileNames = new ArrayList<String>();
		names = new ArrayList<String>();
		countReadsByFile = new ArrayList<CountByKey<Marker>>();
		countBasesByFile = new ArrayList<CountByKey<Marker>>();
		countTypesByFile = new ArrayList<CountByType>();
		markerTypes = new HashMap<Marker, String>();
		coverageByFile = new ArrayList<CoverageByType>();

		if (snpEffectPredictor != null) this.snpEffectPredictor = snpEffectPredictor;
		else this.snpEffectPredictor = new SnpEffectPredictor(new Genome());
	}

	/**
	 * Get marker type
	 * @param marker
	 * @return
	 */
	String markerType(Marker marker) {
		String type = markerTypes.get(marker);
		if (type != null) return type;
		return marker.getClass().getSimpleName();
	}

	/**
	 * Print table to STDOUT
	 */
	public void print() {
		// Show title
		System.out.print("chr\tstart\tend\ttype\tIDs");
		for (int j = 0; j < countReadsByFile.size(); j++)
			System.out.print("\tReads:" + names.get(j) + "\tBases:" + names.get(j));
		System.out.print("\n");

		//---
		// Show counts by marker
		//---
		// Show counts for each marker
		for (Marker key : allMarkers()) {
			// Show 'key' information in first columns
			System.out.print(key.getChromosomeName() //
					+ "\t" + (key.getStart() + 1) //
					+ "\t" + (key.getEnd() + 1) //
					+ "\t" + OutputFormatter.idChain(key) //
			);

			// Show counts for each file
			for (int idx = 0; idx < countReadsByFile.size(); idx++)
				System.out.print("\t" + countReadsByFile.get(idx).get(key) + "\t" + countBasesByFile.get(idx).get(key));
			System.out.print("\n");
		}
		System.out.print("\n");
	}

	/**
	 * Show probabilities
	 * 
	 * @param prob : Probabilities for each 
	 * 
	 * @return A string showing a tab delimited table
	 */
	public String probabilityTable(CountByType prob) {
		StringBuilder sb = new StringBuilder();

		// Create title line
		sb.append("type\tp.binomial"); // Show 'type' information in first columns
		for (int j = 0; j < countReadsByFile.size(); j++)
			sb.append("\treads." + names.get(j) + "\texpected." + names.get(j) + "\tpvalue." + names.get(j));
		sb.append("\n");

		String chrType = Chromosome.class.getSimpleName();

		// Show counts by type
		CountByType countByType = countMarkerTypes(allMarkers());
		for (String type : countByType.keysSorted()) {
			sb.append(type); // Show 'type' information in first columns

			// Binomial probability model
			double p = 0;
			if (prob.contains(type)) {
				p = prob.getScore(type);
				sb.append("\t" + p);
			} else sb.append("\t\t");

			// Show counts for each file
			for (int idx = 0; idx < countReadsByFile.size(); idx++) {
				CountByType count = countTypesByFile.get(idx);

				// Stats
				int n = (int) count.get(chrType); // Number of reads in the file
				int k = (int) count.get(type); // Number of reads hitting this marker type

				long expected = Math.round(count.get(chrType) * p);
				double pvalue = Binomial.get().cdfUpEq(p, k, n);

				if (prob.contains(type)) sb.append("\t" + countTypesByFile.get(idx).get(type) + "\t" + expected + "\t" + pvalue);
				else sb.append("\t" + countTypesByFile.get(idx).get(type) + "\t\t");
			}
			sb.append("\n");
		}

		return sb.toString();
	}

	public void setMarkerTypes(HashMap<Marker, String> markerTypes) {
		this.markerTypes = markerTypes;
	}

	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}

}
