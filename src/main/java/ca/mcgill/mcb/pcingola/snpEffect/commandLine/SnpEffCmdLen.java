package ca.mcgill.mcb.pcingola.snpEffect.commandLine;

import java.util.HashSet;
import java.util.Set;

import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Markers;
import ca.mcgill.mcb.pcingola.probablility.RandMarker;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.SnpEffectPredictor;
import ca.mcgill.mcb.pcingola.stats.CountByType;
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
	CountByType countBases; // Number of bases covered by each marker type
	CountByType countMarkers; // Number of markers (for each marker type)
	CountByType rawCountMarkers; // Number of markers (before join or overlap)
	CountByType rawCountBases; // Number of bases covered by each marker type (befoew join or overlap)
	CountByType prob; // Binomial probability (for each marker type)
	SnpEffectPredictor snpEffectPredictor;

	public SnpEffCmdLen() {
		super();
		countBases = new CountByType();
		countMarkers = new CountByType();
		rawCountMarkers = new CountByType();
		rawCountBases = new CountByType();
		prob = new CountByType();
	}

	/**
	 * Count bases covered for each marker type
	 */
	void countBases() {
		//---
		// Add all markers
		//---
		Markers markers = new Markers();
		markers.add(snpEffectPredictor.getMarkers());
		for (Gene gene : snpEffectPredictor.getGenome().getGenes()) {
			markers.add(gene);
			markers.add(gene.markers());
		}

		// Add all markers (raw counts)
		for (Marker m : markers) {
			String mtype = m.getClass().getSimpleName();
			rawCountMarkers.inc(mtype);
			rawCountBases.inc(mtype, m.size());
		}

		// Count number of bases for each marker type
		for (String mtype : rawCountMarkers.keysSorted()) {
			if (verbose) System.err.print(mtype + ":");

			for (Chromosome chr : snpEffectPredictor.getGenome())
				countBases(mtype, chr, markers);

			if (verbose) System.err.println("");
		}

		// Show chromosomes length
		String mtype = Chromosome.class.getSimpleName();
		for (Chromosome chr : snpEffectPredictor.getGenome()) {
			countBases.inc(mtype, chr.size());
			countMarkers.inc(mtype);
		}
	}

	/**
	 * Count number of bases, for a given chromosome and marker type
	 * @param mtype
	 * @param chr
	 * @param markers
	 * @return
	 */
	void countBases(String mtype, Chromosome chr, Markers markers) {
		String chrName = chr.getChromosomeName();
		if (verbose) System.err.print(" " + chrName);

		// Initialize
		byte busy[] = new byte[chr.size()];
		for (int i = 0; i < busy.length; i++)
			busy[i] = 0;

		for (Marker m : markers) {
			// Same marker type & same chromo? Count bases
			if (m.getChromosomeName().equals(chrName) && m.getClass().getSimpleName().equals(mtype)) {
				for (int i = m.getStart(); i <= m.getEnd(); i++)
					busy[i] = 1;
			}
		}

		int latest = 0;
		for (int i = 0; i < busy.length; i++) {
			// Transition? Count another marker
			if ((i > 0) && (busy[i] != 0) && (busy[i - 1] == 0)) {
				if ((i - latest) <= readLength) countBases.inc(mtype, i - latest); // Intervals are less than one read away? Unify them
				else countMarkers.inc(mtype);
			}

			// Base busy? Count another base
			if (busy[i] != 0) {
				countBases.inc(mtype);
				latest = i;
			}
		}
	}

	public CountByType getCountBases() {
		return countBases;
	}

	public CountByType getCountMarkers() {
		return countMarkers;
	}

	public CountByType getProb() {
		return prob;
	}

	public CountByType getRawCountBases() {
		return rawCountBases;
	}

	public CountByType getRawCountMarkers() {
		return rawCountMarkers;
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
	 * Calculate probabilities
	 */
	void probabilities() {
		// Get total length and count for chromosomes (chromosome size is total genome length)
		String chrType = Chromosome.class.getSimpleName();
		long chrSize = countBases.get(chrType);
		long chrCount = countMarkers.get(chrType);
		if (chrCount <= 0) return; // Zero length genome? Forgot to count bases?

		// Correct readLength 
		int readLength = this.readLength;
		if (readLength < 1) readLength = 1;

		// Probabilities for each marker
		prob = new CountByType();
		for (String mtype : countMarkers.keysSorted()) {
			long size = countBases.get(mtype);
			long count = countMarkers.get(mtype);

			// Calculate and cap probability value
			double p = ((double) (size + (readLength - 1) * count)) / ((double) (chrSize - (readLength - 1) * chrCount));
			p = Math.min(1.0, p);
			p = Math.max(0.0, p);

			prob.addScore(mtype, p);
		}

	}

	/**
	 * Sample and calculate the probability of hitting each type 
	 * of marker (marker.class). Creates 'numReads' reads of 
	 * size 'readLen' and count how many of them hit each marker 
	 * type.
	 */
	CountByType randomSampling(int readLen, int numReads) {
		CountByType countReads = new CountByType();
		RandMarker randMarker = new RandMarker(snpEffectPredictor.getGenome());

		for (int i = 0; i < numReads; i++) {
			// Random read
			Marker read = randMarker.rand(readLen);

			// Where does it hit?
			Set<Marker> regions = snpEffectPredictor.regionsMarkers(read);
			HashSet<String> doneRegion = new HashSet<String>();
			for (Marker m : regions) {
				String clazz = m.getClass().getSimpleName();
				if (!doneRegion.contains(clazz)) {
					countReads.inc(clazz); // Count reads
					doneRegion.add(clazz); // Do not count twice
				}
			}
		}

		return countReads;
	}

	/**
	 * Sample and calculate the probability of hitting each type 
	 * of marker (marker.class). Creates 'numReads' reads of 
	 * size 'readLen' and count how many of them hit each marker 
	 * type. Iterate 'iterations' times to obtain a distribution.
	 */
	void randomSampling(int iterations, int readLen, int numReads) {
		System.out.print("Iteration");
		for (String type : rawCountMarkers.keysSorted())
			System.out.print("\t" + type);
		System.out.println("");

		for (int it = 0; it < iterations; it++) {
			CountByType count = randomSampling(readLen, numReads);
			System.out.print(it);
			for (String type : rawCountMarkers.keysSorted())
				System.out.print("\t" + count.get(type));
			System.out.println("");
		}
	}

	/**
	 * Run
	 * @return 
	 */
	@Override
	public boolean run() {
		if (snpEffectPredictor == null) {
			if (verbose) Timer.showStdErr("Loading config");
			Config config = new Config(genomeVer);

			if (verbose) Timer.showStdErr("Loading predictor");
			snpEffectPredictor = config.loadSnpEffectPredictor();

			if (verbose) Timer.showStdErr("Building interval forest");
			snpEffectPredictor.buildForest();
		}

		if (verbose) Timer.showStdErr("Counting bases");
		countBases(); // Count 
		probabilities(); // Calculate probabilities
		if (!quiet) System.out.println(this);

		// Perform some random sampling
		if ((numIterations > 0) && (readLength > 0)) randomSampling(numIterations, readLength, numReads);

		return true;
	}

	public void setNumIterations(int numIterations) {
		this.numIterations = numIterations;
	}

	public void setNumReads(int numReads) {
		this.numReads = numReads;
	}

	public void setReadLength(int readLength) {
		this.readLength = readLength;
	}

	public void setSnpEffectPredictor(SnpEffectPredictor snpEffectPredictor) {
		this.snpEffectPredictor = snpEffectPredictor;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("marker\tsize\tcount\traw_size\traw_count\tbinomial_p\n");

		probabilities();
		for (String mtype : countMarkers.keysSorted())
			sb.append(mtype + "\t" + countBases.get(mtype) + "\t" + countMarkers.get(mtype) + "\t" + rawCountBases.get(mtype) + "\t" + rawCountMarkers.get(mtype) + "\t" + prob.getScore(mtype) + "\n");

		return sb.toString();
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
