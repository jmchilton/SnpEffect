package ca.mcgill.mcb.pcingola.snpEffect;

import java.util.HashSet;
import java.util.List;

import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Genome;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.interval.SpliceSite;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * Analyze if a set of effects are can create a "Loss Of Function" in a protein.
 * 
 * TODO: What are we supposed to do in cases like UTR_5_DELETED or UTR_3_DELETED?
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
	 * It is assumed that even with a protein coding change at the 
	 * last 5% of the protein, the protein could still be functional.
	 */
	public static final double IGNORE_PROTEIN_CODING_AFTER = 0.95;

	/**
	 *  It is assumed that even with a protein coding change at the 
	 *  first 5% of the protein: 
	 *  	"..suggesting some disrupted transcripts are 
	 *  	rescued by transcriptional reinitiation at an 
	 *  	alternative start codon."
	 */
	public static final double IGNORE_PROTEIN_CODING_BEFORE = 0.05;

	/** 
	 * Larger deletions removing either the first exon or more than 
	 * 50% of the protein-coding sequence of the affected transcript
	 */
	public static final double DELETE_CODING_AFFECTED = 0.50;

	HashSet<Transcript> transcripts;
	HashSet<Gene> genes;

	public LossOfFunction(Genome genome) {
		transcripts = new HashSet<Transcript>();
		genes = new HashSet<Gene>();
	}

	/**
	 * Is this single change a LOF?
	 * @param changeEffect
	 * @return
	 */
	protected boolean isLof(ChangeEffect changeEffect) {
		// Deletion?
		if (changeEffect.getSeqChange().isDel()) return isLofDeletion(changeEffect);

		// The following effect types can be considered LOF
		switch (changeEffect.getEffectType()) {
		case SPLICE_SITE_ACCEPTOR:
		case SPLICE_SITE_DONOR:
			Gpr.debug("SPLICE: " + changeEffect.getMarker());
			// Core splice sites are considered LOF
			if ((changeEffect.getMarker() != null) && (changeEffect.getMarker() instanceof SpliceSite)) {
				// Get splice site marker and check if it is 'core'
				SpliceSite spliceSite = (SpliceSite) changeEffect.getMarker();
				if (spliceSite.isCoreSpliceSite()) return true;

			}
			break;

		case STOP_GAINED:
		case FRAME_SHIFT:
			// It is assumed that even with a protein coding change at the last 5% of the protein, the protein could still be functional.
			double perc = percentCds(changeEffect);
			Gpr.debug("PERCENT: " + perc + "\t" + IGNORE_PROTEIN_CODING_AFTER);
			return (IGNORE_PROTEIN_CODING_BEFORE <= perc) && (perc <= IGNORE_PROTEIN_CODING_AFTER);

		case RARE_AMINO_ACID:
			// This one is not in the referenced papers, but we can assume that RARE AA changes are damaging.
			return true;

		default: // All others are not considered LOF
		}

		return false;
	}

	/**
	 * Can this collection of effects produce a "Loss of function" 
	 * @param changeEffects
	 * @return
	 */
	public boolean isLof(List<ChangeEffect> changeEffects) {
		int lofCount = 0;

		// Iterate over all changeEffects
		for (ChangeEffect changeEffect : changeEffects) {
			if (isLof(changeEffect)) lofCount++;

			transcripts.add(changeEffect.getTranscript()); // Unique transcripts affected (WARNING: null will be added)
			genes.add(changeEffect.getGene()); // Unique genes affected (WARNING: null will be added)
		}

		return lofCount > 0;
	}

	/**
	 * Is this deletion a LOF?
	 * 
	 * Criteria:
	 * 		- More than 50% of coding sequence deleted
	 * 		- First (coding) exon deleted
	 * 
	 * @param changeEffect
	 * @return
	 */
	protected boolean isLofDeletion(ChangeEffect changeEffect) {
		Transcript tr = changeEffect.getTranscript();
		if (tr == null) throw new RuntimeException("Transcript not found for change:\n\t" + changeEffect);

		// Find seqChange
		SeqChange seqChange = changeEffect.getSeqChange();

		// Find coding part of the transcript (i.e. no UTRs)
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
		return (percDeleted > DELETE_CODING_AFFECTED);
	}

	public String lof() {
		return "";
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
	double percentOfTranscriptsAffected(Gene gene) {
		if (gene == null) return 0;

		// Count how many transcript are affected in each gene
		int countAffected = 0;
		for (Transcript tr : gene)
			if (transcripts.contains(tr)) countAffected++;

		return countAffected / ((double) gene.numChilds());
	}
}
