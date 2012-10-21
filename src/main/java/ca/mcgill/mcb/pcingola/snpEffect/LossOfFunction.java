package ca.mcgill.mcb.pcingola.snpEffect;

import java.util.HashSet;
import java.util.List;

import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.interval.SpliceSite;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;
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
	HashSet<Transcript> transcripts;
	HashSet<Gene> genes;

	public LossOfFunction() {
		transcripts = new HashSet<Transcript>();
		genes = new HashSet<Gene>();

		// Config parameters
		config = Config.get();
		ignoreProteinCodingBefore = config.getLofIgnoreProteinCodingBefore();
		ignoreProteinCodingAfter = config.getLofIgnoreProteinCodingAfter();
		deleteProteinCodingBases = config.getLofDeleteProteinCodingBases();
	}

	/**
	 * Is this single change a LOF?
	 * @param changeEffect
	 * @return
	 */
	protected boolean isLof(ChangeEffect changeEffect) {
		// Is this change affecting a protein coding gene?
		Gene gene = changeEffect.getGene();
		if ((gene == null) // No gene affected?
				|| (!gene.isProteinCoding() && !config.isTreatAllAsProteinCoding()) // Not a protein coding gene?
		) return false;

		// Deletion? Is another method to check
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

				// TODO: CHECK IF seqChange hits CORE SITE
				throw new RuntimeException("CHECK IF seqChange hits CORE SITE!");

			}
			break;

		case STOP_GAINED:
		case FRAME_SHIFT:
			// It is assumed that even with a protein coding change at the last 5% of the protein, the protein could still be functional.
			double perc = percentCds(changeEffect);
			Gpr.debug("PERCENT: " + perc + "\t" + ignoreProteinCodingAfter);
			return (ignoreProteinCodingBefore <= perc) && (perc <= ignoreProteinCodingAfter);

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

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();

		for (Gene gene : genes) {
			if (sb.length() > 0) sb.append(','); // Separate by comma
			sb.append(String.format("%s|%s|%d|%.2f", gene.getGeneName(), gene.getId(), gene.numChilds(), percentOfTranscriptsAffected(gene)));
		}

		return sb.toString();
	}
}
