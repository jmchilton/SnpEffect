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
import ca.mcgill.mcb.pcingola.util.Tuple;

/**
 * Calculate the maximum interval length by type, for all markers in a genome
 * 
 * 
 * @author pcingola
 */
public class SnpEffCmdLen extends SnpEff {

	int readLength;
	CountByType markerType;
	SnpEffectPredictor snpEffectPredictor;

	public SnpEffCmdLen() {
		super();
		markerType = new CountByType();
	}

	/**
	 * Count bases ocupied for each marker type
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

		// Add all marker types
		for (Marker m : markers)
			markerType.inc(m.getClass().getSimpleName());
		markerType.remove(Chromosome.class.getSimpleName()); // We don't need chromosomes

		// Count number of bases for each marker type
		System.out.println("count\tsize\ttype");
		for (String mtype : markerType.keysSorted()) {
			if (verbose) System.err.println(mtype);

			long countBases = 0, countMarkers = 0;
			for (Chromosome chr : snpEffectPredictor.getGenome()) {
				Tuple<Long, Long> counters = countBases(mtype, chr, markers);
				countBases += counters.first;
				countMarkers += counters.second;
			}

			System.out.println(countMarkers + "\t" + countBases + "\t" + mtype);
		}

		// Show chromosomes length
		long countBases = 0, countMarkers = 0;
		for (Chromosome chr : snpEffectPredictor.getGenome()) {
			countBases += chr.size();
			countMarkers++;
		}
		System.out.println(countMarkers + "\t" + countBases + "\t" + Chromosome.class.getSimpleName());
	}

	/**
	 * Count number of bases, for a given chromosome and marker type
	 * @param mtype
	 * @param chr
	 * @param markers
	 * @return
	 */
	Tuple<Long, Long> countBases(String mtype, Chromosome chr, Markers markers) {
		long countBases = 0, countMarkers = 0;

		String chrName = chr.getChromosomeName();
		if (verbose) System.err.println("\tChromosome " + chrName);

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
				if ((i - latest) <= readLength) {
					// Intervals are less than one read away? Unify them
					//Gpr.debug("CLOSE: " + mtype + "\t" + (i - latest));
					countBases += (i - latest);
				} else countMarkers++;
			}

			// Base busy? Count another base
			if (busy[i] != 0) {
				countBases++;
				latest = i;
			}

		}

		return new Tuple<Long, Long>(countBases, countMarkers);
	}

	@Override
	public void parseArgs(String[] args) {
		this.args = args;
		for (int i = 0; i < args.length; i++) {

			// Argument starts with '-'?
			if (args[i].equals("-r")) {
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
		snpEffectPredictor = config.loadSnpEffectPredictor();

		if (verbose) Timer.showStdErr("Building interval forest");
		snpEffectPredictor.buildForest();

		if (verbose) Timer.showStdErr("Counting bases");
		countBases();

		if (readLength > 0) {
			// Perform some random sampling
			int numIterations = 10000;
			sample(numIterations, readLength, 1000);
			//			sample(numIterations, readLength, 10000);
			//			sample(numIterations, readLength, 100000);
			//			sample(numIterations, readLength, 1000000);
		}

		return true;
	}

	/**
	 * Sample and calculate the probability of hitting a 'clazz' marker when  'numReads' reads of size 'readLen'
	 * @param readLen
	 */
	CountByType sample(int readLen, int numReads) {
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
	 * Sampling random reads.
	 * 
	 * @param countByType
	 */
	void sample(int iterations, int readLen, int numReads) {
		System.out.print("Iteration");
		for (String type : markerType.keysSorted())
			System.out.print("\t" + type);
		System.out.println("");

		for (int it = 0; it < iterations; it++) {
			CountByType count = sample(readLen, numReads);
			System.out.print(it);
			for (String type : markerType.keysSorted())
				System.out.print("\t" + count.get(type));
			System.out.println("");
		}

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
		System.err.println("\t-r <num> : Assume a read size of 'num' bases.");
		System.exit(-1);
	}

}
