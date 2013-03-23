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

	SnpEffectPredictor snpEffectPredictor;
	HashSet<Marker> done;
	CountByType count, size;
	Markers maxMarkers;
	int readLength;

	public SnpEffCmdLen() {
		super();
		done = new HashSet<Marker>();
		count = new CountByType();
		size = new CountByType();
		maxMarkers = new Markers();
	}

	/**
	 * Add this marker to the final list of markers
	 * @param maxMarker
	 */
	void add(Marker maxMarker) {
		maxMarkers.add(maxMarker); // Add

		String clazz = maxMarker.getClass().getSimpleName(); // Some stats
		count.inc(clazz);
		size.inc(clazz, maxMarker.size());

		// Show
		if (verbose) System.out.println(clazz + "\t" + maxMarker.getChromosomeName() + "\t" + (maxMarker.getStart() + 1) + "\t" + (maxMarker.getEnd() + 1));
	}

	/**
	 * Expand marker by 'readLength' bases on each size
	 * @param m
	 * @param readLength
	 * @return
	 */
	Marker expand(Marker m, int readLength) {
		if (readLength == 0) return m;
		return expand(m, m.getStart(), m.getEnd(), readLength);
	}

	/**
	 * Expand marker by 'readLength' bases on each size
	 * @param m
	 * @param readLength
	 * @return
	 */
	Marker expand(Marker m, int start, int end, int readLength) {
		Marker mexp = m.clone();

		start = Math.max(start - readLength, 0);
		mexp.setStart(start);

		end = end + readLength;
		mexp.setEnd(end);

		return mexp;
	}

	/**
	 * Calculate length for all markers
	 * @param m
	 */
	void len() {
		// For all Genes and sub intervals
		for (Chromosome chr : snpEffectPredictor.getGenome())
			len(chr);

		// For all Genes and sub intervals
		for (Gene gene : snpEffectPredictor.getGenome().getGenes()) {
			len(gene);
			for (Marker m : gene.markers())
				len(m);
		}

		// All other intervals
		for (Marker m : snpEffectPredictor.getMarkers())
			len(m);

	}

	/**
	 * Calculate length for a marker
	 * @param m
	 */
	void len(Marker m) {
		if (!done.contains(m)) {
			if (debug) System.out.println(m.toStr());
			Marker maxm = maxIntersect(m);
			if (debug && (m.size() != maxm.size())) System.out.println("\t" + m.size() + "\t" + maxm.size() + "\t" + maxm.toStr());

			add(m); // Add to set
			done.add(m); // Make sure we don't use this marker again
		}
	}

	/**
	 * Get a marker representing the maximum intercept of this marker
	 * @param mori
	 * @return
	 */
	Marker maxIntersect(Marker mori) {
		Marker mexp = expand(mori, readLength); // Expand by readLength

		Markers markers = snpEffectPredictor.intersects(mexp); // Get markers intersecting expanded marker

		// Iterate until there is no size change 
		Marker mprev = mori, mnew = null;
		while (true) {
			mnew = maxIntersect(mprev, mexp, markers);
			if (mnew.size() == mprev.size()) return mnew;
			if (debug) Gpr.debug("\tPrevious size: " + mprev.size() + "\tnew size: " + mnew.size());
			mprev = mnew;
		}
	}

	/**
	 * Get a marker representing the maximum of Markers (same class a 'm')
	 * @param m
	 * @return
	 */
	Marker maxIntersect(Marker mori, Marker mexp, Markers markers) {
		int start = mori.getStart(), end = mori.getEnd();
		boolean update = false;

		for (Marker m : markers) {
			// Intersects a gene?
			if (m instanceof Gene) {
				// Look for all sub-markers
				Markers markersGene = ((Gene) m).markers();

				for (Marker mm : markersGene) {
					// Same class? Intersects expanded marker? => Add
					if ((mm.getClass() == mori.getClass()) && (mexp.intersects(mm))) {
						if (debug) Gpr.debug("\tMerge:" + mori.toStr() + "\t" + mm.toStr());
						start = Math.min(start, mm.getStart());
						end = Math.max(end, mm.getEnd());
						update = true;
						mexp = expand(mori, start, end, readLength); // Update expanded marker
						done.add(mm);// Make sure we don't use the same marker again
					}
				}
			} else if (m.getClass() == mori.getClass()) {
				if (debug) Gpr.debug("\tMerge:" + mori.toStr() + "\t" + m.toStr());
				start = Math.min(start, m.getStart());
				end = Math.max(end, m.getEnd());
				update = true;
				mexp = expand(mori, start, end, readLength); // Update expanded marker
				done.add(m);// Make sure we don't use the same marker again
			}
		}

		// Updated? Create new marker using new coordinates
		if (update) {
			Marker mnew = mori.clone();
			mnew.setStart(start);
			mnew.setEnd(end);
			return mnew;
		}

		return mori; // No update, return original marker
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
		Timer.showStdErr("Loading config");
		Config config = new Config(genomeVer);

		Timer.showStdErr("Loading predictor");
		snpEffectPredictor = config.loadSnpEffectPredictor();

		Timer.showStdErr("Building interval forest");
		snpEffectPredictor.buildForest();

		Timer.showStdErr("Calculating (max intersection) interval length");
		len();

		// Show totals
		System.out.println("\nTotals:\n\tcount\tsize\ttype");
		for (String clazz : count.keySet())
			System.out.println("\t" + count.get(clazz) + "\t" + size.get(clazz) + "\t" + clazz);

		// Perform some random sampling
		//sample(count, 10000, 100, 1000);
		// sample(count, 10000, 100, 10000);
		// sample(count, 10000, 100, 100000);

		return true;
	}

	/**
	 * Sampling random reads.
	 * 
	 * @param countByType
	 */
	void sample(CountByType countByType, int iterations, int readLen, int numReads) {
		System.out.print("Iteration");
		for (String type : countByType.keysSorted())
			System.out.print("\t" + type);
		System.out.println("");

		for (int it = 0; it < iterations; it++) {
			CountByType count = sample(readLen, numReads);
			System.out.print(it);
			for (String type : countByType.keysSorted())
				System.out.print("\t" + count.get(type));
			System.out.println("");
		}

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
			HashSet<String> done = new HashSet<String>();
			for (Marker m : regions) {
				String clazz = m.getClass().getSimpleName();
				if (!done.contains(clazz)) {
					countReads.inc(clazz); // Count reads
					done.add(clazz); // Do not count twice
				}
			}
		}

		return countReads;
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
