package ca.mcgill.mcb.pcingola.interval;

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;

import ca.mcgill.mcb.pcingola.interval.codonChange.CodonChange;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.stats.ObservedOverExpectedCpG;

/**
 * Codon position
 * @author pcingola
 */
class CodonPosition {

	public int codonNum = -1;
	public int codonIndex = -1;
}

/**
 * Interval for a transcript, as well as some other information: exons, utrs, cds, etc.
 * 
 * @author pcingola
 */
public class Transcript extends IntervalAndSubIntervals<Exon> {

	private static final long serialVersionUID = -2665025617916107311L;

	ArrayList<Utr> utrs;
	ArrayList<Cds> cdss; // Sometimes we have additional CDS information
	Upstream upstream;
	Downstream downstream;
	String cds = null; // Coding sequence
	String bioType; // Transcript biotype
	int cdsStart = -1;
	int cdsEnd = -1;
	boolean proteinCoding = false; // Is this a protein-coding transcript?

	public Transcript(Marker gene, int start, int end, int strand, String id) {
		super(gene, start, end, strand, id);
		utrs = new ArrayList<Utr>();
		cdss = new ArrayList<Cds>();
		type = EffectType.TRANSCRIPT.toString();
	}

	/**
	 * Add a CDS
	 * @param cdsInt
	 */
	public void add(Cds cdsInt) {
		cdss.add(cdsInt);
		cds = null;
	}

	/**
	 * Add a UTR
	 * @param utr
	 */
	public void add(Utr utr) {
		utrs.add(utr);
		cds = null;
	}

	/**
	 * Add missing UTRs. See utrFromCds() method.
	 * @param missingUtrs
	 */
	boolean addMissingUtrs(Markers missingUtrs, boolean verbose) {
		missingUtrs.sort(false, strand < 0);

		// Get min/max CDS positions
		int minCds = Integer.MAX_VALUE;
		int maxCds = 0;
		for (Cds c : cdss) {
			minCds = Math.min(minCds, c.getStart());
			maxCds = Math.max(maxCds, c.getEnd());
		}

		if (verbose) System.out.println("Transcript '" + id + "' has missing UTRs. Strand: " + strand + " (minCds: " + minCds + " , maxCds: " + maxCds + "):");

		// Add intervals
		boolean retVal = false;
		for (Marker mu : missingUtrs) {
			Exon eint = intersectingExon(mu);
			if (eint == null) throw new RuntimeException("Cannot find exon for UTR: " + mu);
			Utr toAdd = null;

			if (isStrandPlus()) {
				if (mu.getEnd() <= minCds) toAdd = new Utr5prime(eint, mu.getStart(), mu.getEnd(), strand, mu.getId());
				else if (mu.getStart() >= maxCds) toAdd = new Utr3prime(this, mu.getStart(), mu.getEnd(), strand, mu.getId());
			} else {
				if (mu.getStart() >= maxCds) toAdd = new Utr5prime(eint, mu.getStart(), mu.getEnd(), strand, mu.getId());
				else if (mu.getEnd() <= minCds) toAdd = new Utr3prime(this, mu.getStart(), mu.getEnd(), strand, mu.getId());
			}

			// OK?
			if (toAdd != null) {
				add(toAdd);
				if (verbose) System.out.println("\tAdding " + toAdd);
				retVal = true;
			}
		}

		return retVal;
	}

	public boolean adjust() {
		boolean changed = false;
		int strandSumTr = 0;
		int newStart = start, newEnd = end;
		if (newStart == 0 && newEnd == 0) {
			newStart = Integer.MAX_VALUE;
			newEnd = Integer.MIN_VALUE;
		}

		for (Exon exon : sortedStrand()) {
			newStart = Math.min(newStart, exon.getStart());
			newEnd = Math.max(newEnd, exon.getEnd());
			strandSumTr += exon.getStrand(); // Some exons have incorrect strands, we use the strand indicated by most exons
		}

		for (Utr utr : getUtrs()) {
			newStart = Math.min(newStart, utr.getStart());
			newEnd = Math.max(newEnd, utr.getEnd());
		}

		// Change transcript strand?
		int newStrand = strandSumTr >= 0 ? 1 : -1;
		if (strand != newStrand) {
			strand = newStrand;
			changed = true;
		}

		// Change start?
		if (start != newStart) {
			start = newStart;
			changed = true;
		}

		// Change end?
		if (end != newEnd) {
			end = newEnd;
			changed = true;
		}

		return changed;
	}

	/**
	 * Calculate CDS start and CDS end
	 */
	synchronized void calcCdsStartEnd() {
		if (cdsStart < 0) {
			// Calculate coding start (after 5 prime UTR)

			if (utrs.isEmpty()) { // No UTRs => use exons
				cdsStart = (isStrandPlus() ? end : start); // cdsStart is the position of the first base in the CDS (i.e. the first base after all 5'UTR)
				cdsEnd = (isStrandPlus() ? start : end); // cdsEnd is the position of the last base in the CDS (i.e. the first base before all 3'UTR)

				for (Exon ex : this) {
					if (isStrandPlus()) {
						cdsStart = Math.min(cdsStart, ex.getStart());
						cdsEnd = Math.max(cdsEnd, ex.getEnd());
					} else {
						cdsStart = Math.max(cdsStart, ex.getEnd());
						cdsEnd = Math.min(cdsEnd, ex.getStart());
					}
				}

			} else {
				cdsStart = (isStrandPlus() ? start : end); // cdsStart is the position of the first base in the CDS (i.e. the first base after all 5'UTR)
				cdsEnd = (isStrandPlus() ? end : start); // cdsEnd is the position of the last base in the CDS (i.e. the first base before all 3'UTR)

				for (Utr utr : utrs) {
					if (utr instanceof Utr5prime) {
						if (isStrandPlus()) cdsStart = Math.max(cdsStart, utr.getEnd() + 1);
						else cdsStart = Math.min(cdsStart, utr.getStart() - 1);
					} else if (utr instanceof Utr3prime) {
						if (isStrandPlus()) cdsEnd = Math.min(cdsEnd, utr.getStart() - 1);
						else cdsEnd = Math.max(cdsEnd, utr.getEnd() + 1);
					}
				}
			}
		}
	}

	/**
	 * Retrieve coding sequence
	 */
	public String cds() {
		if (cds != null) return cds;

		// Concatenate all exons
		List<Exon> exons = sortedStrand();
		StringBuilder sequence = new StringBuilder();
		int utr5len = 0, utr3len = 0;

		// 5 prime UTR length
		for (Utr utr : get5primeUtrs())
			utr5len += utr.size();

		// Append all exon sequences
		for (Exon eint : exons)
			sequence.append(eint.getSequence());

		// 3 prime UTR length
		for (Utr utr : get3primeUtrs())
			utr3len += utr.size();

		// Cut 5 prime UTR and 3 prime UTR points
		int subEnd = sequence.length() - utr3len;

		if (utr5len > subEnd) cds = "";
		else cds = sequence.substring(utr5len, subEnd);

		return cds;
	}

	/**
	 * Calculate base number in a CDS where 'pos' maps
	 * 
	 * @returns Base number or '-1' if it does not map to a coding base
	 */
	public int cdsBaseNumber(int pos, boolean usePrevBaseIntron) {
		// Doesn't hit this transcript?
		if (!intersects(pos)) return -1;

		// Is it in UTR instead of CDS? 
		if (isUtr(pos)) return -1;

		// Calculate cdsStart and cdsEnd (if not already done)
		calcCdsStartEnd();

		// All exons..
		int firstCdsBaseInExon = 0; // Where the exon maps to the CDS (i.e. which CDS base number does the first base in this exon maps to).
		for (Exon eint : sortedStrand()) {
			if (eint.intersects(pos)) {
				int cdsBaseInExon; // cdsBaseInExon: base number relative to the beginning of the coding part of this exon (i.e. excluding 5'UTRs)
				if (strand >= 0) cdsBaseInExon = pos - Math.max(eint.getStart(), cdsStart);
				else cdsBaseInExon = Math.min(eint.getEnd(), cdsStart) - pos;

				cdsBaseInExon = Math.max(0, cdsBaseInExon);

				return firstCdsBaseInExon + cdsBaseInExon;
			} else {
				// Before exon begins?
				if ((isStrandPlus() && (pos < eint.getStart())) // Before exon begins (positive strand)?
						|| (isStrandMinus() && (pos > eint.getEnd()))) // Before exon begins (negative strand)?
					return firstCdsBaseInExon - (usePrevBaseIntron ? 1 : 0);
			}

			if (isStrandPlus()) firstCdsBaseInExon += Math.max(0, eint.getEnd() - Math.max(eint.getStart(), cdsStart) + 1);
			else firstCdsBaseInExon += Math.max(0, Math.min(cdsStart, eint.getEnd()) - eint.getStart() + 1);
		}

		return firstCdsBaseInExon - 1;
	}

	/**
	 * Calculate chromosome position as function of CDS number 
	 * 
	 * @returns An array mapping 'pos[cdsBaseNumber] = chromosmalPos' 
	 */
	public int[] cdsBaseNumber2ChrPos() {
		calcCdsStartEnd();

		int cds2pos[] = new int[cds().length()];
		for (int i = 0; i < cds2pos.length; i++)
			cds2pos[i] = -1;

		int cdsMin = Math.min(cdsStart, cdsEnd);
		int cdsMax = Math.max(cdsStart, cdsEnd);

		// For each exon, add CDS position to array
		int cdsBaseNum = 0;
		for (Exon exon : sortedStrand()) {
			int min = isStrandPlus() ? exon.getStart() : exon.getEnd();
			int step = isStrandPlus() ? 1 : -1;
			for (int pos = min; exon.intersects(pos); pos += step)
				if ((cdsMin <= pos) && (pos <= cdsMax)) cds2pos[cdsBaseNum++] = pos;
		}

		return cds2pos;
	}

	/**
	 * Analyze SNPs in this transcript.
	 * Add changeEffect to 'changeEffect'
	 */
	public String codonByCdsBaseNumber(int cdsBaseNumber) {
		int codonNum = cdsBaseNumber / CodonChange.CODON_SIZE;
		int min = codonNum * CodonChange.CODON_SIZE;
		int max = codonNum * CodonChange.CODON_SIZE + CodonChange.CODON_SIZE;
		if ((min >= 0) && (max <= cds().length())) return cds().substring(min, max).toUpperCase();
		return null;
	}

	/**
	 * Calculate CpG bias: number of CpG / expected[CpG]
	 * @return
	 */
	public double cpgExonBias() {
		ObservedOverExpectedCpG oe = new ObservedOverExpectedCpG();
		return oe.oe(this);
	}

	/**
	 * Count total CpG in this transcript's exons
	 * @return
	 */
	public int cpgExons() {
		ObservedOverExpectedCpG oe = new ObservedOverExpectedCpG();
		return oe.observed(this);
	}

	/**
	 * Create splice sites for this transcript
	 * @return
	 */
	public List<SpliceSite> createSpliceSites() {
		List<SpliceSite> list = new LinkedList<SpliceSite>();

		// For each gene, transcript and exon
		ArrayList<Exon> exons = new ArrayList<Exon>();
		exons.addAll(sortedStrand());

		if (exons.size() > 0) {
			for (int i = 0; i < exons.size(); i++) {

				Exon exon = exons.get(i);
				Exon prev = (i >= 1 ? exons.get(i - 1) : null);
				Exon next = (i < exons.size() - 1 ? exons.get(i + 1) : null);

				//---
				// Distance to previous exon
				//---
				if (prev != null) {
					int dist = 0;
					if (strand >= 0) dist = exon.getStart() - prev.getEnd() - 1;
					else dist = prev.getStart() - exon.getEnd() - 1;

					// Acceptor splice site: before exon start, but not before first exon
					SpliceSite ss = exon.createSpliceSiteAcceptor(dist);
					if (ss != null) list.add(ss);
				}

				//---
				// Distance to next exon
				//---
				if (next != null) {
					int dist = 0;
					if (strand >= 0) dist = next.getStart() - exon.getEnd() - 1;
					else dist = exon.getStart() - next.getEnd() - 1;

					// Donor splice site: after exon end, but not after last exon
					SpliceSite ss = exon.createSpliceSiteDonor(dist);
					if (ss != null) list.add(ss);
				}

				// Sanity check
				int rank = i + 1;
				if (exon.getRank() != rank) throw new RuntimeException("Rank numbers do not march: " + rank + " != " + exon.getRank());
			}
		}

		return list;
	}

	/**
	 * Creates a list of UP/DOWN stream regions (for each transcript)
	 * Upstream (downstream) stream is defined as upDownLength before (after) transcript
	 */
	public void createUpDownStream(int upDownLength) {
		Chromosome chr = getChromosome();
		int min = chr.getStart(), max = chr.getEnd();

		// Create up/down stream intervals and add them to the list
		if (isStrandPlus()) {
			upstream = new Upstream(this, Math.max(start - upDownLength, min), Math.max(start - 1, min), 1, id);
			downstream = new Downstream(this, Math.min(end + 1, max), Math.min(end + upDownLength, max), 1, id);
		} else {
			upstream = new Upstream(this, Math.min(end + 1, max), Math.min(end + upDownLength, max), 1, id);
			downstream = new Downstream(this, Math.max(start - upDownLength, min), Math.max(start - 1, min), 1, id);
		}
	}

	/**
	 * Deletes redundant exons (i.e. exons that are totally included in other exons).
	 * Does the same for CDSs.
	   Does the same for UTRs.
	 */
	public void deleteRedundant() {
		// Delete redundant exons
		Collection<Marker> toDelete = MarkerUtil.redundant(subintervals());
		for (Marker exon : toDelete)
			subIntervals.remove(exon.getId());

		// Delete redundant CDS
		toDelete = MarkerUtil.redundant(cdss);
		for (Marker cds : toDelete)
			cdss.remove(cds);

		// Delete redundant CDS
		toDelete = MarkerUtil.redundant(utrs);
		for (Marker utr : toDelete)
			utrs.remove(utr);
	}

	/**
	 * Create a list of 3 prime UTRs
	 */
	public List<Utr3prime> get3primeUtrs() {
		ArrayList<Utr3prime> list = new ArrayList<Utr3prime>();
		for (Utr utr : utrs)
			if (utr instanceof Utr3prime) list.add((Utr3prime) utr);
		return list;
	}

	/**
	 * Create a list of 5 prime UTRs
	 */
	public List<Utr5prime> get5primeUtrs() {
		ArrayList<Utr5prime> list = new ArrayList<Utr5prime>();
		for (Utr utr : utrs)
			if (utr instanceof Utr5prime) list.add((Utr5prime) utr);
		return list;
	}

	public String getBioType() {
		return bioType;
	}

	/**
	 * Get all CDSs
	 * @return
	 */
	public List<Cds> getCds() {
		return cdss;
	}

	public int getCdsEnd() {
		calcCdsStartEnd();
		return cdsEnd;
	}

	public int getCdsStart() {
		calcCdsStartEnd();
		return cdsStart;
	}

	public Downstream getDownstream() {
		return downstream;
	}

	public Upstream getUpstream() {
		return upstream;
	}

	/**
	 * Get all UTRs
	 * @return
	 */
	public List<Utr> getUtrs() {
		return utrs;
	}

	/**
	 * Return the first exon that intersects 'interval' (null if not found)
	 * @param interval
	 * @return
	 */
	public Exon intersectingExon(Marker interval) {
		for (Exon ei : this)
			if (ei.intersects(interval)) return ei;
		return null;
	}

	@Override
	protected boolean isAdjustIfParentDoesNotInclude(Marker parent) {
		return true;
	}

	/**
	 * Is this seqChange in the CDS part of this transcript?
	 * @param seqChange
	 * @return
	 */
	boolean isCds(SeqChange seqChange) {
		calcCdsStartEnd();

		int cs = cdsStart;
		int ce = cdsEnd;

		if (isStrandMinus()) {
			cs = cdsEnd;
			ce = cdsStart;
		}

		return (seqChange.getEnd() >= cs) && (seqChange.getStart() <= ce);
	}

	public boolean isProteinCoding() {
		return proteinCoding;
	}

	/**
	 * Does this 'pos' hit a UTR?
	 * @param pos
	 * @return
	 */
	boolean isUtr(int pos) {
		// Is it in UTR instead of CDS? 
		for (Utr utr : utrs)
			if (utr.intersects(pos)) return true;
		return false;
	}

	/**
	 * Retrieve coding sequence AND the UTRs (mRNA = 5'UTR + CDS + 3'UTR)
	 * I.e. Concatenate all exon sequences
	 */
	public String mRna() {
		List<Exon> exons = sortedStrand();

		// Concatenate all exons
		StringBuilder sequence = new StringBuilder();
		for (Exon eint : exons)
			sequence.append(eint.getSequence());

		return sequence.toString();
	}

	public String protein() {
		if (!Config.get().isTreatAllAsProteinCoding() && !isProteinCoding()) return "";
		return codonTable().aa(cds());
	}

	/**
	 * Assign ranks to exons
	 */
	public boolean rankExons() {
		boolean changed = false;
		int rank = 1;
		for (Exon exon : sortedStrand()) {
			if (rank != exon.getRank()) {
				exon.setRank(rank);
				changed = true;
			}
			rank++;
		}
		return changed;
	}

	/**
	 * Get some details about the effect on this transcript
	 * @param seqChange
	 * @return
	 */
	@Override
	public List<ChangeEffect> seqChangeEffect(SeqChange seqChange, ChangeEffect changeEffect) {
		if (!intersects(seqChange)) return ChangeEffect.emptyResults(); // Sanity check

		// Create a list of changes
		ArrayList<ChangeEffect> changeEffectList = new ArrayList<ChangeEffect>();

		//---
		// Hits a UTR region?
		//---
		boolean includedInUtr = false;
		for (Utr utr : utrs)
			if (utr.intersects(seqChange)) {
				// Calculate the effect
				List<ChangeEffect> chEffList = utr.seqChangeEffect(seqChange, changeEffect.clone());
				if (!chEffList.isEmpty()) changeEffectList.addAll(chEffList);

				// Is this seqChange fully included in the UTR?
				includedInUtr |= utr.includes(seqChange);
			}

		// Since the effect was fully included in the UTR, we are done.
		if (includedInUtr) return changeEffectList;

		//---
		// Analyze non-coding transcripts (or 'interval' seqChanges)
		//---
		if ((!Config.get().isTreatAllAsProteinCoding() && !isProteinCoding()) || seqChange.isInterval()) {
			// Do we have exon information for this transcript?
			if (!subintervals().isEmpty()) {
				// Add all exons
				for (Exon exon : this) {
					if (exon.intersects(seqChange)) {
						ChangeEffect cheff = changeEffect.clone();
						cheff.set(exon, EffectType.EXON, "");
						changeEffectList.add(cheff);
					}
				}

				// Did not hit an exon or a UTR? => Must hit an intron
				if (changeEffectList.isEmpty()) {
					ChangeEffect cheff = changeEffect.clone();
					cheff.set(this, EffectType.INTRON, "");
					changeEffectList.add(cheff);
				}
			} else {
				// No exons annotated? Just mark it as hitting a transcript
				ChangeEffect cheff = changeEffect.clone();
				cheff.set(this, EffectType.TRANSCRIPT, "");
				changeEffectList.add(cheff);
			}

			return changeEffectList;
		}

		//---
		// This is a protein coding transcript.
		// We analyze codon replacement effect
		//---
		if (isCds(seqChange)) {
			// Get codon change effect 
			CodonChange codonChange = new CodonChange(seqChange, this, changeEffect);
			List<ChangeEffect> codonChangesList = codonChange.calculate();

			// Get exons. In most cases it's only one exon per change.
			for (ChangeEffect resCodon : codonChangesList) {
				for (Exon exon : this) {
					if (exon.intersects(seqChange)) {
						ChangeEffect chEff = resCodon.clone();
						changeEffectList.addAll(exon.seqChangeEffect(seqChange, chEff));
					}
				}
			}
		}

		// No exons? => It's intronic region
		if (changeEffectList.isEmpty()) {
			changeEffect.set(this, EffectType.INTRON, "");
			return changeEffect.newList();
		}

		return changeEffectList;
	}

	public void setBioType(String bioType) {
		this.bioType = bioType;
	}

	public void setProteinCoding(boolean proteinCoding) {
		this.proteinCoding = proteinCoding;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();

		sb.append(getChromosomeName() + ":" + start + "-" + end);
		sb.append(", strand: " + (strand >= 0 ? "+" : "-"));
		if ((id != null) && (id.length() > 0)) sb.append(", id:" + id);
		if (isProteinCoding()) sb.append(", Protein");

		if (numChilds() > 0) {
			sb.append("\n");
			for (Utr utr : get5primeUtrs())
				sb.append("\t\t5'UTR   :\t" + utr + "\n");

			sb.append("\t\tExons:\n");
			for (Exon eint : sorted())
				sb.append("\t\t" + eint + "\n");

			for (Utr utr : get3primeUtrs())
				sb.append("\t\t3'UTR   :\t" + utr + "\n");

			// We may show CDS
			if (isProteinCoding()) {
				sb.append("\t\tCDS     :\t" + cds() + "\n");
				sb.append("\t\tProtein :\t" + protein() + "\n");
			}
		}

		return sb.toString();
	}

	/**
	 * Calculate UTR regions from CDSs
	 */
	public boolean utrFromCds(boolean verbose) {
		if (cdss.size() <= 0) return false; // Cannot do this if we don't have CDS information

		// All exons minus all UTRs and CDS should give us the missing UTRs
		Markers exons = new Markers();
		Markers minus = new Markers();

		// Add all exons
		for (Exon e : this)
			exons.add(e);

		// Add all UTRs and CDSs to the 'minus' set
		for (Utr uint : getUtrs())
			minus.add(uint);

		for (Cds cint : cdss)
			minus.add(cint);

		Markers missingUtrs = exons.minus(minus); // Perform interval minus
		if (missingUtrs.size() > 0) return addMissingUtrs(missingUtrs, verbose); // Anything left? => There was a missing UTR
		return false;
	}

}
