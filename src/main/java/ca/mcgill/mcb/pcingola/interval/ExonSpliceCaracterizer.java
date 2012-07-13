package ca.mcgill.mcb.pcingola.interval;

import java.util.HashMap;

import ca.mcgill.mcb.pcingola.interval.Exon.ExonSpliceType;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.SnpEffectPredictor;
import ca.mcgill.mcb.pcingola.stats.CountByType;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * Caracterize exons based on alternative splicing
 * 
 * References: "Alternative splicing and evolution - diversification, exon definition and function"  (see Box 1)
 * 
 * @author pablocingolani
 */
public class ExonSpliceCaracterizer {

	boolean verbose = false;
	Genome genome;
	HashMap<Exon, Exon.ExonSpliceType> typeByExon;
	CountByType countByType = new CountByType();

	public ExonSpliceCaracterizer(Genome genome) {
		this.genome = genome;
	}

	public ExonSpliceCaracterizer(String genomeVer) {
		if (verbose) Timer.showStdErr("Loading dataabse");
		Config config = new Config(genomeVer);
		SnpEffectPredictor snpEffectPredictor = config.loadSnpEffectPredictor();
		genome = snpEffectPredictor.getGenome();
		if (verbose) Timer.showStdErr("done.");
	}

	/**
	 * Caracterize all exons
	 */
	public CountByType caracterize(boolean verbose) {
		this.verbose = verbose;

		if (verbose) Timer.showStdErr("Run");
		typeByExon = new HashMap<Exon, Exon.ExonSpliceType>();
		type();
		if (verbose) Timer.showStdErr("Done.\n\tTotal exons  : " + countExons() + "\n\tExons marked : " + typeByExon.size() + "\n" + countByType);
		return countByType;
	}

	/**
	 * Count number of exons
	 * @return
	 */
	int countExons() {
		int count = 0;
		for (Gene g : genome.getGenes())
			for (Transcript tr : g)
				count += tr.numChilds();
		return count;
	}

	/**
	 * Does the marker intersect any exon in 'tr'?
	 * @param m
	 * @param tr
	 * @return
	 */
	boolean intersectsAnyExon(Marker m, Transcript tr) {
		for (Exon e : tr)
			if (m.intersects(e)) return true;
		return false;
	}

	/**
	 * Is thins an ALTERNATIVE_3SS exon?
	 * @param exon
	 * @param gene
	 * @return
	 */
	boolean isAlt3ss(Exon exon, Gene gene) {
		for (Transcript tr : gene)
			for (Exon e : tr) {
				if (exon.intersects(e)) {
					if (exon.isStrandPlus()) {
						// Same exon end, different exon start?
						if ((exon.getStart() != e.getStart()) && (exon.getEnd() == e.getEnd())) return true;
					} else {
						// Same exon end, different exon start? (negative strand)
						if ((exon.getStart() == e.getStart()) && (exon.getEnd() != e.getEnd())) return true;
					}
				}
			}

		return false;
	}

	/**
	 * Is thins an ALTERNATIVE_5SS exon?
	 * @param exon
	 * @param gene
	 * @return
	 */
	boolean isAlt5ss(Exon exon, Gene gene) {
		for (Transcript tr : gene)
			for (Exon e : tr) {
				if (exon.intersects(e)) {
					if (exon.isStrandPlus()) {
						// Same exon start, different exon end?
						if ((exon.getStart() == e.getStart()) && (exon.getEnd() != e.getEnd())) return true;
					} else {
						// Same exon start, different exon end? (negative strand)
						if ((exon.getStart() != e.getStart()) && (exon.getEnd() == e.getEnd())) return true;
					}
				}
			}

		return false;
	}

	/**
	 * Is this exon mutually exclusive with another exon?
	 * @param exon
	 * @param gene
	 * @return
	 */
	boolean isMutEx(Exon exon, Gene gene) {
		if (gene.numChilds() <= 1) return false;

		//---
		// Make a list of all unique exons 
		//---
		String exonKey = key(exon);
		HashMap<String, Exon> uniqEx = new HashMap<String, Exon>();
		for (Transcript tr : gene)
			for (Exon e : tr) {
				String ekey = key(e);
				if (!exonKey.equals(ekey)) uniqEx.put(ekey, e);
			}

		//---
		// For each unque exon, compare if it is mutually exclusive with 'exon'
		//---
		Transcript exonTr = (Transcript) exon.getParent();
		for (Exon e : uniqEx.values()) {
			ExonSpliceType type = typeByExon.get(e);;
			if (type == Exon.ExonSpliceType.SKIPPED) { // Only check for these exons (avoid all ALT*)
				boolean xor = true;
				for (Transcript tr : gene) {
					if (exonTr.intersects(tr) && exon.intersects(tr) && e.intersects(exonTr)) { // Make sure both exons intersect both transcripts (otherwise they cannot be mutually exclusive)
						xor &= intersectsAnyExon(e, tr) ^ intersectsAnyExon(exon, tr);
					} else xor = false;
				}

				// XOR is true? => Mutually exclusive
				if (xor) return true;
			}
		}

		return false;
	}

	/**
	 * Create a simple hash key based on choromosomal position
	 * @param m
	 * @return
	 */
	String key(Marker m) {
		return m.getChromosomeName() + ":" + m.getStart() + "-" + m.getEnd();
	}

	/**
	 * Mark exons types
	 */
	void type() {
		// Find retained exons
		for (Gene g : genome.getGenes()) {
			// Count exons
			CountByType count = new CountByType();
			for (Transcript tr : g)
				for (Exon e : tr)
					count.inc(key(e));

			// Label exons
			int countTr = g.numChilds();
			for (Transcript tr : g) {
				for (Exon e : tr) {
					String eKey = key(e);
					int countEx = (int) count.get(eKey);

					// Is this exon maintained in all transcripts? 
					if (countEx == countTr) type(e, Exon.ExonSpliceType.RETAINED);
					else {
						if (isAlt3ss(e, g)) type(e, Exon.ExonSpliceType.ALTTENATIVE_3SS);
						else if (isAlt5ss(e, g)) type(e, Exon.ExonSpliceType.ALTTENATIVE_5SS);
						else if (tr.numChilds() > 1) {
							if (e.getRank() == 1) type(e, Exon.ExonSpliceType.ALTTENATIVE_PROMOMOTER);
							else if (e.getRank() == tr.numChilds()) {
								Gpr.debug("ALTTENATIVE_POLY_A: " + e.getId() + "\n" + tr);
								type(e, Exon.ExonSpliceType.ALTTENATIVE_POLY_A);
							} else type(e, Exon.ExonSpliceType.SKIPPED);
						}
					}
				}
			}
		}

		// Now analyze if there are mutually exclusive exons
		for (Gene g : genome.getGenes())
			for (Transcript tr : g) {
				for (Exon e : tr) {
					ExonSpliceType type = typeByExon.get(e);
					if (type == ExonSpliceType.SKIPPED) { // Try to re-annotate only these
						if (isMutEx(e, g)) type(e, Exon.ExonSpliceType.MUTUALLY_EXCLUSIVE);
					}
				}
			}

		Timer.showStdErr("done.");
	}

	/**
	 * Mark this exons as 'type'
	 * @param e
	 * @param type
	 */
	void type(Exon e, Exon.ExonSpliceType type) {
		e.spliceType = type;
		countByType.inc(type.toString());
		typeByExon.put(e, type);
	}
}
