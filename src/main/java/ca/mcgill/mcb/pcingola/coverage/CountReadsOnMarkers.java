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
import ca.mcgill.mcb.pcingola.stats.plot.GoogleBarChart;
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
	List<String> fileNames;
	List<String> names;
	ArrayList<CountByKey<Marker>> countReadsByFile;
	ArrayList<CountByKey<Marker>> countBasesByFile;
	ArrayList<CountByType> countTypesByFile;
	SnpEffectPredictor snpEffectPredictor;
	MarkerTypes markerTypes;
	ArrayList<CoverageByType> coverageByFile;
	CountByType readsByFile;
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
		fileNames.add(samFileName);
		names.add(Gpr.removeExt(Gpr.baseName(samFileName)));
	}

	public void addMarkerType(Marker marker, String type) {
		markerTypes.addType(marker, type);
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
		for (String fileName : fileNames) {
			try {
				if (verbose) Timer.showStdErr("Reading file '" + fileName + "'");
				CountByKey<Marker> countReads = new CountByKey<Marker>();
				CountByKey<Marker> countBases = new CountByKey<Marker>();
				CountByType countTypes = new CountByType();
				CoverageByType coverageByType = new CoverageByType();

				countFile(fileName, countReads, countBases, countTypes, coverageByType);

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
				Timer.showStdErr("Finished reding file " + fileName + "\n\tTotal reads: " + readsByFile.get(fileName));
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
				countMarker(fileName, read, countReads, countBases, countTypes, coverageByType);
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
	void countMarker(String fileName, Marker read, CountByKey<Marker> countReads, CountByKey<Marker> countBases, CountByType countTypes, CoverageByType coverageByType) {
		// Find all intersects
		Markers regions = snpEffectPredictor.queryDeep(read);

		// Count total reads
		readsByFile.inc(fileName);

		// Count each marker
		HashSet<String> doneClass = new HashSet<String>();
		for (Marker m : regions) {
			countReads.inc(m); // Count reads
			countBases.inc(m, m.intersectSize(read)); // Count number bases that intersect

			// Count by marker type (make sure we only count once per read)
			String type = markerTypes.getType(m);
			String subtype = markerTypes.getSubType(m);

			if (!doneClass.contains(type)) {
				countTypes.inc(type); // Count reads
				doneClass.add(type); // Do not count twice

				PosStats posStats = coverageByType.getOrCreate(type);
				posStats.sample(read, m);
			}

			// Count sub-type if any
			if ((subtype != null) && !doneClass.contains(subtype)) {
				countTypes.inc(subtype); // Count reads
				doneClass.add(subtype); // Do not count twice

				PosStats posStats = coverageByType.getOrCreate(subtype);
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
		for (Marker marker : markersToCount) {
			String type = markerTypes.getType(marker);
			String subtype = markerTypes.getSubType(marker);
			countByMarkerType.inc(type);
			if (subtype != null) countByMarkerType.inc(subtype);
		}
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

						countMarker(fileName, read, countReads, countBases, countTypes, coverageByType);
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
				countMarker(fileName, read, countReads, countBases, countTypes, coverageByType);
				if (verbose) Gpr.showMark(readNum, SHOW_EVERY);
				readNum++;
			} catch (Exception e) {
				countExceptions++;
				if (countExceptions < 10) e.printStackTrace();
				else if (countExceptions == 10) System.err.println("Not showing more exceptions!");
			}
		}
	}

	public MarkerTypes getMarkerTypes() {
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

	/**
	 * Show charts in html 
	 * @return
	 */
	public String html() {
		StringBuilder sbHead = new StringBuilder();
		StringBuilder sbBody = new StringBuilder();

		//---
		// Barchart: By Marker types (all files toghether)
		//---

		// Create 3 charts: One for all intervals, one for Exons and one for Introns.
		HashSet<String> keySetAll = new HashSet<String>();
		HashSet<String> keySetExon = new HashSet<String>();
		HashSet<String> keySetIntron = new HashSet<String>();
		for (CountByType ct : countTypesByFile)
			for (String key : ct.keySet())
				if (key.startsWith("Exon:")) keySetExon.add(key);
				else if (key.startsWith("Intron:")) keySetIntron.add(key);
				else keySetAll.add(key);

		keySetAll.remove("Chromosome"); // We don't want this number in the chart (usually it's too big)
		HashMap<String, HashSet<String>> keySets = new HashMap<String, HashSet<String>>();
		keySets.put("", keySetAll);
		keySets.put("Exons", keySetExon);
		keySets.put("Introns", keySetIntron);

		// Sort names
		ArrayList<String> keySetNames = new ArrayList<String>();
		keySetNames.addAll(keySets.keySet());
		Collections.sort(keySetNames);

		// Create one barchart for each keySet
		for (String ksname : keySetNames) {
			HashSet<String> keySet = keySets.get(ksname);
			// Sort keys
			ArrayList<String> keys = new ArrayList<String>();
			keys.addAll(keySet);
			Collections.sort(keys);

			// Add all columns
			GoogleBarChart barchart = new GoogleBarChart("Count by file " + ksname);
			barchart.setxLables(keys);

			GoogleBarChart barchartPercent = new GoogleBarChart("Count by file " + ksname + " [Percent]");
			barchartPercent.setxLables(keys);

			// Add all files
			for (int i = 0; i < names.size(); i++) {
				String name = names.get(i);
				CountByType ct = countTypesByFile.get(i);

				// Add all values
				ArrayList<String> columnValues = new ArrayList<String>();
				for (String key : keys)
					if (keys.contains(key)) columnValues.add(ct.get(key) + "");

				// Add column to chart
				barchart.addColumn(name, columnValues);
				barchartPercent.addColumn(name, columnValues);
			}

			// Add header and body
			sbHead.append(barchart.toStringHtmlHeader());
			sbBody.append(barchart.toStringHtmlBody());
			barchartPercent.percentColumns();
			sbHead.append(barchartPercent.toStringHtmlHeader());
			sbBody.append(barchartPercent.toStringHtmlBody());
		}

		//---
		// Barchart: By Marker types (one by file)
		//---

		// Add all files
		for (int i = 0; i < names.size(); i++) {
			String name = names.get(i);

			// Create one barchart for each keySet
			for (String ksname : keySetNames) {
				// Sort keys
				HashSet<String> keySet = keySets.get(ksname);
				ArrayList<String> keys = new ArrayList<String>();
				keys.addAll(keySet);
				Collections.sort(keys);

				GoogleBarChart barchart = new GoogleBarChart("Count by file " + name + " " + ksname);
				barchart.setxLables(keys);

				// Add all columns
				CountByType ct = countTypesByFile.get(i);

				// Add all values
				ArrayList<String> columnValues = new ArrayList<String>();
				for (String key : keys)
					if (keys.contains(key)) columnValues.add(ct.get(key) + "");

				// Create chart
				barchart.addColumn(name, columnValues);
				sbHead.append(barchart.toStringHtmlHeader());
				sbBody.append(barchart.toStringHtmlBody());
			}
		}

		//---
		// Genomic region charts
		//---
		ArrayList<GoogleGeneRegionChart> genRegCharts = new ArrayList<GoogleGeneRegionChart>();
		for (int i = 0; i < names.size(); i++) {
			String name = names.get(i);
			CoverageByType cvt = coverageByFile.get(i);

			GoogleGeneRegionChart grc = new GoogleGeneRegionChart(cvt, name);
			genRegCharts.add(grc);
		}

		// Add all headers
		for (GoogleGeneRegionChart grc : genRegCharts)
			sbHead.append(grc.toStringHtmlHeader());

		// Add all bodies
		for (GoogleGeneRegionChart grc : genRegCharts)
			sbBody.append(grc.toStringHtmlBody());

		// Return all html code
		return "<head>\n" + sbHead.toString() + "\n</head>\n" + sbBody.toString();
	}

	/**
	 * Initialize
	 * @param snpEffectPredictor
	 */
	void init(SnpEffectPredictor snpEffectPredictor) {
		fileNames = new ArrayList<String>();
		names = new ArrayList<String>();
		countReadsByFile = new ArrayList<CountByKey<Marker>>();
		countBasesByFile = new ArrayList<CountByKey<Marker>>();
		countTypesByFile = new ArrayList<CountByType>();
		markerTypes = new MarkerTypes();
		coverageByFile = new ArrayList<CoverageByType>();
		readsByFile = new CountByType();

		if (snpEffectPredictor != null) this.snpEffectPredictor = snpEffectPredictor;
		else this.snpEffectPredictor = new SnpEffectPredictor(new Genome());
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
			if ((prob != null) && prob.contains(type)) {
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

				if ((prob != null) && prob.contains(type)) sb.append("\t" + countTypesByFile.get(idx).get(type) + "\t" + expected + "\t" + pvalue);
				else sb.append("\t" + countTypesByFile.get(idx).get(type) + "\t\t");
			}
			sb.append("\n");
		}

		return sb.toString();
	}

	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}

	/**
	 * Print table to STDOUT
	 */
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();

		// Show title
		sb.append("chr\tstart\tend\ttype\tIDs");
		for (int j = 0; j < countReadsByFile.size(); j++)
			sb.append("\tReads:" + names.get(j) + "\tBases:" + names.get(j));
		sb.append("\n");

		//---
		// Show counts by marker
		//---
		// Show counts for each marker
		for (Marker key : allMarkers()) {
			// Show 'key' information in first columns
			sb.append(key.getChromosomeName() //
					+ "\t" + (key.getStart() + 1) //
					+ "\t" + (key.getEnd() + 1) //
					+ "\t" + OutputFormatter.idChain(key) //
			);

			// Show counts for each file
			for (int idx = 0; idx < countReadsByFile.size(); idx++)
				sb.append("\t" + countReadsByFile.get(idx).get(key) + "\t" + countBasesByFile.get(idx).get(key));
			sb.append("\n");
		}
		sb.append("\n");

		return sb.toString();
	}

}
