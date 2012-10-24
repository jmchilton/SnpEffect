package ca.mcgill.mcb.pcingola.snpEffect;

import java.util.Collection;
import java.util.HashSet;

import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.interval.SpliceSite;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;

/**
 * Analyze if a set of effects are can create a "Loss Of Function" 
 * and "Nonsense mediated decays" effects.
 * 
 * TODO: Add branch points? (We have to analyze correlation with expression)
 * 
 * TODO: Other NextProt markers? (We have to analyze correlation with expression)
 * 
 * TODO: What are we supposed to do in cases like UTR_5_DELETED or UTR_3_DELETED?
 *       So far we are considering them as moderate impact. 
 * 
 * Of course, this is a prediction based on analysis 
 * of groups of "putative effects". Proper wet-lab 
 * validation is required to infer "real" LOF. 
 * 
 * References: I used the LOF definition used in the 
 * following paper "A Systematic Survey of Loss-of-Function 
 * Variants in Human Protein-Coding Genes", Science, 2012
 * 
 * From the paper: 
 * 		We adopted a definition for LoF variants
 * 		expected to correlate with complete loss of function 
 * 		of the affected transcripts: stop codon-introducing 
 * 		(nonsense) or splice site-disrupting single-nucleotide 
 * 		variants (SNVs), insertion/deletion (indel) variants 
 * 		predicted to disrupt a transcript's reading frame, or 
 * 		larger deletions removing either the first exon or more 
 * 		than 50% of the protein-coding sequence of the affected 
 * 		transcript.
 * 
 * 		Both nonsense SNVs and frameshift indels are enriched toward the 3' end
 * 		of the affected gene, consistent with a greater tolerance to truncation 
 * 		close to the end of the coding sequence (Fig. 1C); putative LoF variants
 * 		identified in the last 5% of the coding region were thus systematically 
 * 		removed from our high-confidence set.
 * 
 * 
 * @author pcingola
 */
public class LossOfFunction {

	/**
	 * Number of bases before last exon-exon junction that nonsense 
	 * mediated decay is supposed to occur 
	 */
	public static final int MND_BASES_BEFORE_LAST_JUNCTION = 50;

	/** 
	 * It is assumed that even with a protein coding change at the 
	 * last 5% of the protein, the protein could still be functional.
	 */
	public static final double DEFAULT_IGNORE_PROTEIN_CODING_AFTER = 0.95;
	public double ignoreProteinCodingAfter;

	/**
	 *  It is assumed that even with a protein coding change at the 
	 *  first 5% of the protein: 
	 *  	"..suggesting some disrupted transcripts are 
	 *  	rescued by transcriptional reinitiation at an 
	 *  	alternative start codon."
	 */
	public static final double DEFAULT_IGNORE_PROTEIN_CODING_BEFORE = 0.05;
	public double ignoreProteinCodingBefore;

	/** 
	 * Larger deletions removing either the first exon or more than 
	 * 50% of the protein-coding sequence of the affected transcript
	 */
	public static final double DEFAULT_DELETE_PROTEIN_CODING_BASES = 0.50;
	public double deleteProteinCodingBases;

	Config config;
	HashSet<Transcript> transcriptsLof;
	HashSet<Gene> genesLof;
	HashSet<Transcript> transcriptsNmd;
	HashSet<Gene> genesNmd;
	Collection<ChangeEffect> changeEffects;
	int lofCount = -1; // Number of loss of function effects
	int nmdCount = -1; // Number of nonsense mediated decay effects

	public LossOfFunction(Collection<ChangeEffect> changeEffects) {
		this.changeEffects = changeEffects;
		transcriptsLof = new HashSet<Transcript>();
		genesLof = new HashSet<Gene>();
		transcriptsNmd = new HashSet<Transcript>();
		genesNmd = new HashSet<Gene>();

		// Config parameters
		config = Config.get();
		ignoreProteinCodingBefore = config.getLofIgnoreProteinCodingBefore();
		ignoreProteinCodingAfter = config.getLofIgnoreProteinCodingAfter();
		deleteProteinCodingBases = config.getLofDeleteProteinCodingBases();
	}

	/**
	 * Can this collection of effects produce a "Loss of function" 
	 * @param changeEffects
	 * @return
	 */
	public boolean isLof() {
		// Need to calculate?
		if (lofCount < 0) {
			lofCount = 0;

			// Iterate over all changeEffects
			for (ChangeEffect changeEffect : changeEffects)
				if (isLof(changeEffect)) lofCount++;
		}

		return lofCount > 0;
	}

	/**
	 * Is this single change a LOF?
	 * 
	 * Criteria:
	 * 		1) Core splice sites acceptors or donors (only CORE ones)
	 * 		2) Stop gained (if this happens at the last part of the protein, we assume it has no effect)
	 * 		3) Frame shifts
	 * 
	 * @param changeEffect
	 * @return
	 */
	protected boolean isLof(ChangeEffect changeEffect) {
		// Is this change affecting a protein coding gene?
		Gene gene = changeEffect.getGene();
		if ((gene == null) // No gene affected?
				|| (!gene.isProteinCoding() && !config.isTreatAllAsProteinCoding()) // Not a protein coding gene?
		) return false;

		// Is this variant a LOF?
		boolean lof = false;
		if (changeEffect.getSeqChange().isDel()) lof = isLofDeletion(changeEffect); // Deletion? Is another method to check
		else {
			// The following effect types can be considered LOF
			switch (changeEffect.getEffectType()) {
			case SPLICE_SITE_ACCEPTOR:
			case SPLICE_SITE_DONOR:
				// Core splice sites are considered LOF
				if ((changeEffect.getMarker() != null) && (changeEffect.getMarker() instanceof SpliceSite)) {
					// Get splice site marker and check if it is 'core'
					SpliceSite spliceSite = (SpliceSite) changeEffect.getMarker();
					if (spliceSite.intersectsCoreSpliceSite(changeEffect.getSeqChange())) lof = true; // Does it intersect the CORE splice site?
				}
				break;

			case STOP_GAINED:
				lof = isNmd(changeEffect);
				break;

			case FRAME_SHIFT:
				// It is assumed that even with a protein coding change at the last 5% of the protein, the protein could still be functional.
				double perc = percentCds(changeEffect);
				lof = (ignoreProteinCodingBefore <= perc) && (perc <= ignoreProteinCodingAfter);
				break;

			case RARE_AMINO_ACID:
			case START_LOST:
				// This one is not in the referenced papers, but we assume that RARE AA and START_LOSS changes are damaging.
				lof = true;
				break;

			default: // All others are not considered LOF
				lof = false;
			}
		}

		// Update sets
		if (lof) {
			transcriptsLof.add(changeEffect.getTranscript()); // Unique transcripts affected (WARNING: null will be added)
			genesLof.add(changeEffect.getGene()); // Unique genes affected (WARNING: null will be added)
		}

		return lof;
	}

	/**
	 * Is this deletion a LOF?
	 * 
	 * Criteria:
	 * 		1) First (coding) exon deleted
	 * 		2) More than 50% of coding sequence deleted
	 * 
	 * @param changeEffect
	 * @return
	 */
	protected boolean isLofDeletion(ChangeEffect changeEffect) {
		Transcript tr = changeEffect.getTranscript();
		if (tr == null) throw new RuntimeException("Transcript not found for change:\n\t" + changeEffect);

		//---
		// Criteria:
		// 		1) First (coding) exon deleted
		//---
		if (changeEffect.getEffectType() == EffectType.EXON_DELETED) {
			Exon exon = changeEffect.getExon();
			if (exon == null) throw new RuntimeException("Cannot retrieve 'exon' from EXON_DELETED effect!");
			if (tr.getFirstCodingExon() == exon) return true;
		}

		//---
		// Criteria:
		// 		2) More than 50% of coding sequence deleted
		//---

		// Find coding part of the transcript (i.e. no UTRs)
		SeqChange seqChange = changeEffect.getSeqChange();
		int cdsStart = tr.getCdsStart();
		int cdsEnd = tr.getCdsEnd();
		Marker coding = new Marker(seqChange.getChromosome(), cdsStart, cdsEnd, 1, "");

		// Create an interval intersecting the CDS and the deletion
		int start = Math.max(cdsStart, seqChange.getStart());
		int end = Math.min(cdsEnd, seqChange.getEnd());
		if (start >= end) return false; // No intersections with coding part of the exon? => not LOF
		Marker codingDeleted = new Marker(seqChange.getChromosome(), start, end, 1, "");

		// Count:
		//   - number of coding bases deleted 
		//   - number of coding bases
		int codingBasesDeleted = 0, codingBases = 0;
		for (Exon exon : tr) {
			codingBasesDeleted += codingDeleted.intersectSize(exon);
			codingBases += coding.intersectSize(exon);
		}

		// More than a threshold? => It is a LOF
		double percDeleted = codingBasesDeleted / ((double) codingBases);
		return (percDeleted > deleteProteinCodingBases);
	}

	/**
	 * Can this collection of effects produce a "Nonsense mediated decay"?
	 * @param changeEffects
	 * @return
	 */
	public boolean isNmd() {
		if (nmdCount < 0) isLof(); // Need to calculate?
		return nmdCount > 0;
	}

	/**
	 * Is this single change a LOF?
	 * 
	 * Criteria:
	 * 		1) Core splice sites acceptors or donors (only CORE ones)
	 * 		2) Stop gained (if this happens at the last part of the protein, we assume it has no effect)
	 * 		3) Frame shifts
	 * 
	 * @param changeEffect
	 * @return
	 */
	protected boolean isNmd(ChangeEffect changeEffect) {
		Transcript tr = changeEffect.getTranscript();
		if (tr == null) throw new RuntimeException("Transcript not found for change:\n\t" + changeEffect);

		// Only one exon? Nothing to do (there is no exon-exon junction)
		if (tr.numChilds() <= 1) return false;

		// Find last valid NMD position
		int lastNmdPos = lastNmdPos(tr);

		// Does this change affect the region 'before' this last NMD position? => It is assumed to be NMD
		SeqChange seqChange = changeEffect.getSeqChange();

		boolean nmd;
		if (tr.isStrandPlus()) nmd = seqChange.getStart() <= lastNmdPos;
		else nmd = lastNmdPos <= seqChange.getEnd();

		// Update sets and counters
		if (nmd) {
			transcriptsNmd.add(changeEffect.getTranscript()); // Unique transcripts affected (WARNING: null will be added)
			genesNmd.add(changeEffect.getGene()); // Unique genes affected (WARNING: null will be added)
			nmdCount++;
		}

		return nmd;
	}

	/**
	 * Find the last position where a nonsense mediated decay is supposed to occurr
	 * This is 50 bases (MND_BASES_BEFORE_LAST_JUNCTION bases) before the last exon-exon junction.
	 * 
	 * @param tr
	 * @return
	 */
	protected int lastNmdPos(Transcript tr) {
		//---
		// Get last exon
		//---
		int cdsEnd = tr.getCdsEnd();
		Exon lastExon = null;
		for (Exon exon : tr.sortedStrand())
			if (exon.intersects(cdsEnd)) lastExon = exon;

		// Sanity check
		if (lastExon == null) throw new RuntimeException("Cannot find last coding exon for transcript '" + tr.getId() + "' (cdsEnd: " + cdsEnd + ")\n\t" + tr);

		//---
		// Find that position of MND_BASES_BEFORE_LAST_JUNCTION before the last exon-exon junction
		//---
		int lastExonJunction = tr.isStrandPlus() ? lastExon.getStart() : lastExon.getEnd();
		int chrPos[] = tr.cdsBaseNumber2ChrPos();
		int lastNmdPos = -1;
		for (int cdsi = chrPos.length - 1; cdsi >= 0; cdsi--) {
			if (chrPos[cdsi] == lastExonJunction) {
				if (cdsi >= MND_BASES_BEFORE_LAST_JUNCTION) lastNmdPos = chrPos[cdsi - MND_BASES_BEFORE_LAST_JUNCTION - 1];
				else return 0; // Out of CDS range
				return lastNmdPos;
			}
		}

		throw new RuntimeException("Cannot find last exon junction position for transcript '" + tr.getId() + "'\n\t" + tr);
		// return -1;
	}

	/**
	 * Which percentile of the protein does this effect hit?
	 * @param changeEffect
	 * @return
	 */
	double percentCds(ChangeEffect changeEffect) {
		int cdsLen = changeEffect.getAaLength();
		int codonNum = changeEffect.getCodonNum();
		if ((cdsLen >= 0) && (codonNum >= 0)) return codonNum / ((double) cdsLen);
		return Double.NaN;
	}

	/**
	 * What percentile of the transcripts in this gene are affected?
	 * @param gene
	 * @return 
	 */
	double percentOfTranscriptsAffected(Gene gene, HashSet<Transcript> transcripts) {
		if (gene == null) return 0;

		// Count how many transcript are affected in each gene
		int countAffected = 0;
		for (Transcript tr : gene)
			if (transcripts.contains(tr)) countAffected++;

		return countAffected / ((double) gene.numChilds());
	}

	@Override
	public String toString() {
		return (isLof() ? "LOF=" + vcfLofValue() + " " : "") //
				+ (isNmd() ? "NMD=" + vcfNmdValue() : "") //
		;
	}

	/**
	 * Get LOF value for VCF info field
	 * @return
	 */
	public String vcfLofValue() {
		StringBuilder sb = new StringBuilder();

		for (Gene gene : genesLof) {
			if (sb.length() > 0) sb.append(','); // Separate by comma
			sb.append(String.format("%s|%s|%d|%.2f", gene.getGeneName(), gene.getId(), gene.numChilds(), percentOfTranscriptsAffected(gene, transcriptsLof)));
		}

		return sb.toString();
	}

	/**
	 * Get NMD value for VCF info field
	 * @return
	 */
	public String vcfNmdValue() {
		StringBuilder sb = new StringBuilder();

		for (Gene gene : genesNmd) {
			if (sb.length() > 0) sb.append(','); // Separate by comma
			sb.append(String.format("%s|%s|%d|%.2f", gene.getGeneName(), gene.getId(), gene.numChilds(), percentOfTranscriptsAffected(gene, transcriptsNmd)));
		}

		return sb.toString();
	}

}