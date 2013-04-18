package ca.mcgill.mcb.pcingola.snpEffect.commandLine;

import ca.mcgill.mcb.pcingola.fileIterator.BedFileIterator;
import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Intergenic;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Markers;
import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.interval.Utr;
import ca.mcgill.mcb.pcingola.outputFormatter.OutputFormatter;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.SnpEffectPredictor;
import ca.mcgill.mcb.pcingola.util.Timer;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;

/**
 * Command line: Find closes marker to each variant
 * 
 * @author pcingola
 */
public class SnpEffCmdClosest extends SnpEff {

	public static final String CLOSEST = "CLOSEST";
	public static final String INFO_LINE = "##INFO=<ID=" + CLOSEST + ",Number=4,Type=String,Description=\"Closest exon: Distance (bases), exons Id, transcript Id, gene name\">";

	boolean bedFormat = false;
	String inFile = "";
	SnpEffectPredictor snpEffectPredictor;

	//	IntervalForest intervalForest;

	public SnpEffCmdClosest() {
		super();
		command = "closestExon";
	}

	public SnpEffCmdClosest(Config config) {
		super();
		command = "closestExon";
		this.config = config;
		inFile = config.getFileNameProteins();
	}

	/**
	 * Update header
	 * @param vcf
	 */
	void addHeaderLines(VcfFileIterator vcf) {
		vcf.getVcfHeader().addLine("##SnpEffVersion=\"" + SnpEff.VERSION + "\"");
		vcf.getVcfHeader().addLine("##SnpEffCmd=\"" + commandLineStr(false) + "\"");
		vcf.getVcfHeader().addLine(INFO_LINE);
	}

	/**
	 * Iterate over VCF file, find closest exons and annotate vcf lines
	 */
	void bedIterate() {
		// Open file
		BedFileIterator bfi = new BedFileIterator(inFile, config.getGenome(), 0);
		bfi.setCreateChromos(true); // Any 'new' chromosome in the input file will be created (otherwise an error will be thrown)

		for (SeqChange bed : bfi) {
			try {
				// Find closest exon
				Marker closestMarker = findClosestMarker(bed);

				String id = bed.getId();

				// Update INFO fields if any exon was found
				if (closestMarker != null) {
					int dist = closestMarker.distance(bed);
					id = (id.isEmpty() ? "" : bed.getId() + ";") + dist + "," + OutputFormatter.idChain(closestMarker, ",", false);
				}

				// Show output
				System.out.println(bed.getChromosomeName() //
						+ "\t" + bed.getStart() // BED format: Zero-based position
						+ "\t" + (bed.getEnd() + 1) // BED format: End base is not included
						+ "\t" + id //
				);

			} catch (Exception e) {
				e.printStackTrace(); // Show exception and move on...
			}
		}
	}

	/**
	 * Find closest marker
	 * @param queryMarker
	 */
	Marker findClosestMarker(Marker queryMarker) {
		int initialExtension = 1000;

		Chromosome chr = queryMarker.getChromosome();
		if ((chr != null) && (chr.size() > 0)) {

			// Extend interval to capture 'close' markers
			for (int extend = initialExtension; extend < chr.size(); extend *= 2) {
				int start = Math.max(queryMarker.getStart() - extend, 0);
				int end = queryMarker.getEnd() + extend;
				Marker extended = new Marker(chr, start, end, 1, "");

				// Find all markers that intersect with 'extended interval'
				Markers markers = snpEffectPredictor.queryDeep(extended);
				Marker closest = findClosestMarker(queryMarker, markers);
				if (closest != null) return closest;
			}
		}

		// Nothing found
		return null;
	}

	/**
	 * Find closest marker to query (in markers collection)
	 * @param queryMarker
	 * @param markers
	 * @return
	 */
	Marker findClosestMarker(Marker queryMarker, Markers markers) {
		// We prefer the closest marker from the longest transcript)
		int minDist = Integer.MAX_VALUE;
		int maxTrLen = 0;

		Marker minDistMarker = null;
		for (Marker m : markers) {
			// We don't care about these
			if ((m instanceof Chromosome) || (m instanceof Intergenic) || (m instanceof Gene) || (m instanceof Transcript)) continue;

			// Find closest marker
			int dist = m.distance(queryMarker);
			if (dist <= minDist) {

				Transcript tr = findTranscript(m);
				if (tr != null) {
					// Find closest marker in largest transcript
					int trLen = tr.mRna().length();
					if (trLen > maxTrLen) {
						maxTrLen = trLen;
						minDist = dist;
						minDistMarker = m;
					} else if (trLen == maxTrLen) {
						// Prefer Utr to Exon (more descriptive)
						if ((minDistMarker instanceof Exon) && (m instanceof Utr)) minDistMarker = m;
					}
				}
			}
		}

		return minDistMarker;
	}

	Transcript findTranscript(Marker m) {
		if (m instanceof Transcript) return (Transcript) m;
		return (Transcript) m.findParent(Transcript.class);
	}

	/**
	 * Parse command line arguments
	 */
	@Override
	public void parseArgs(String[] args) {
		this.args = args;
		for (int i = 0; i < args.length; i++) {

			// Argument starts with '-'?
			if (args[i].startsWith("-")) {
				if (args[i].equals("-bed")) bedFormat = true;
				else usage("Unknow option '" + args[i] + "'");
			} else if (genomeVer.isEmpty()) genomeVer = args[i];
			else if (inFile.isEmpty()) inFile = args[i];
			else usage("Unknow parameter '" + args[i] + "'");
		}

		// Check: Do we have all required parameters?
		if (genomeVer.isEmpty()) usage("Missing genomer_version parameter");
		if (inFile.isEmpty()) usage("Missing protein_file parameter");
	}

	/**
	 * Run command
	 */
	@Override
	public boolean run() {
		// Load config
		if (config == null) {
			if (verbose) Timer.showStdErr("Reading configuration...");
			config = new Config(genomeVer, configFile); // Read configuration
			if (verbose) Timer.showStdErr("done");
		}

		if (verbose) Timer.showStdErr("Loading predictor...");
		config.loadSnpEffectPredictor();
		if (verbose) Timer.showStdErr("done");

		if (verbose) Timer.showStdErr("Building interval forest...");
		snpEffectPredictor = config.getSnpEffectPredictor();
		snpEffectPredictor.buildForest();
		if (verbose) Timer.showStdErr("done");

		if (verbose) Timer.showStdErr("Reading file '" + inFile + "'");
		if (bedFormat) bedIterate();
		else vcfIterate();
		if (verbose) Timer.showStdErr("done");

		return true;
	}

	@Override
	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}

	/**
	 * Show usage and exit
	 */
	@Override
	public void usage(String message) {
		if (message != null) System.err.println("Error: " + message + "\n");
		System.err.println("snpEff version " + SnpEff.VERSION);
		System.err.println("Usage: snpEff closestExon [options] genome_version file.vcf");
		System.err.println("\nOptions:");
		System.err.println("\t-bed          : Input format is BED. Default: VCF");
		System.exit(-1);
	}

	/**
	 * Iterate over VCF file, find closest exons and annotate vcf lines
	 */
	void vcfIterate() {
		// Open file
		VcfFileIterator vcf = new VcfFileIterator(inFile, config.getGenome());
		vcf.setCreateChromos(true); // Any 'new' chromosome in the input file will be created (otherwise an error will be thrown)

		boolean header = true;
		for (VcfEntry ve : vcf) {
			try {
				if (header) {
					// Update and show header
					addHeaderLines(vcf);
					String headerStr = vcf.getVcfHeader().toString();
					if (!headerStr.isEmpty()) System.out.println(headerStr);
					header = false;
				}

				// Find closest exon
				Marker closestMarker = findClosestMarker(ve);

				// Update INFO fields if any exon was found
				if (closestMarker != null) {
					int dist = closestMarker.distance(ve);
					String value = dist + "," + OutputFormatter.idChain(closestMarker, ",", false);
					ve.addInfo(CLOSEST, value);
				}

				// Show output
				System.out.println(ve);
			} catch (Exception e) {
				e.printStackTrace(); // Show exception and move on...
			}
		}
	}

}
