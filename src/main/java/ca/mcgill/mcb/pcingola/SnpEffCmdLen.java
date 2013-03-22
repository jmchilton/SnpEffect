package ca.mcgill.mcb.pcingola;

import java.util.HashSet;

import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Markers;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.SnpEffectPredictor;
import ca.mcgill.mcb.pcingola.snpEffect.commandLine.SnpEff;
import ca.mcgill.mcb.pcingola.stats.CountByType;
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
	 * Make a complete markers list, by adding all gene's child markers 
	 * @param markers
	 * @return
	 */
	Markers complete(Marker marker, Markers markers) {
		Markers newMarkers = new Markers();
		for (Marker m : markers) {

			// Is this a gene?
			if (m instanceof Gene) {
				// Look for all sub-markers
				Markers markersGene = ((Gene) m).markers();

				for (Marker mm : markersGene) {
					// Same class? Intersect? => Add
					if ((mm.getClass() == marker.getClass()) && (marker.intersects(mm))) {
						newMarkers.add(mm);
						done.add(mm);// Make sure we don't use this marker again
					}
				}
			} else if (m.getClass() == marker.getClass()) {
				newMarkers.add(m); // Add all markers of the same class
				done.add(m); // Make sure we don't use this marker again
			}
		}
		return newMarkers;
	}

	/**
	 * Get a marker representing the maximum intercept of this marker
	 * @param m
	 * @return
	 */
	Marker max(Marker m) {
		Markers markers = snpEffectPredictor.intersects(m); // Get markers intersecting this one
		if (markers.size() <= 1) return m; // Nothing else?
		markers = complete(m, markers); // Add all gene's sub-markers
		return max(m, markers); // Get max spanning interval
	}

	/**
	 * Get a marker representing the maximum of Markers (same class a 'm')
	 * @param m
	 * @return
	 */
	Marker max(Marker m, Markers markers) {
		int start = m.getStart(), end = m.getEnd();

		for (Marker mm : markers) {
			if (mm.getStart() < start) {
				start = mm.getStart();
				if (debug) System.out.println("\tMERGE:\t" + m.toStr() + "\t" + mm.toStr());
			}

			if (end < mm.getEnd()) {
				end = mm.getEnd();
				if (debug) System.out.println("\tMERGE:\t" + m.toStr() + "\t" + mm.toStr());
			}
		}

		return new Marker(m.getParent(), start, end, 1, m.getId());
	}

	@Override
	public void parseArgs(String[] args) {
		this.args = args;
		for (int i = 0; i < args.length; i++) {

			// Argument starts with '-'?
			if (genomeVer.isEmpty()) genomeVer = args[i];
			else usage("Unknow parameter '" + args[i] + "'");
		}

		// Check: Do we have all required parameters?
		if (genomeVer.isEmpty()) usage("Missing genomer_version parameter");
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

		Timer.showStdErr("Calculating interval sizes");

		// For all Genes and sub intervals
		for (Gene gene : snpEffectPredictor.getGenome().getGenes()) {
			for (Marker m : gene.markers()) {
				if (!done.contains(m)) {
					if (debug) System.out.println(m.toStr());
					Marker maxm = max(m);
					if (debug && (m.size() != maxm.size())) System.out.println("\t" + m.size() + "\t" + maxm.size() + "\t" + maxm.toStr());

					done.add(m); // Make sure we don't use this marker again
					add(m);
				}
			}
		}

		// All other intervals
		for (Marker m : snpEffectPredictor.getMarkers()) {
			if (!done.contains(m)) {
				if (debug) System.out.println(m.toStr());
				Marker maxm = max(m);
				if (debug && (m.size() != maxm.size())) System.out.println("\t" + m.size() + "\t" + maxm.size() + "\t" + maxm.toStr());

				done.add(m); // Make sure we don't use this marker again
				add(m);
			}
		}

		// Show totals
		System.out.println("\nTotals:\n\tcount\tsize\ttype");
		for (String clazz : count.keySet())
			System.out.println("\t" + count.get(clazz) + "\t" + size.get(clazz) + "\t" + clazz);

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
		System.exit(-1);
	}

}
