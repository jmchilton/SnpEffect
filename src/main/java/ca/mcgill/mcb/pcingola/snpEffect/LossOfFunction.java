package ca.mcgill.mcb.pcingola.snpEffect;

import java.util.List;

import ca.mcgill.mcb.pcingola.interval.Genome;
import ca.mcgill.mcb.pcingola.interval.SpliceSite;

/**
 * Analyze if a set of effects are can create a 
 * "Loss Of Function" in a protein.
 * 
 * Of course, this is a prediction based on analysis 
 * of groups of "putative effects". Proper wet-lab 
 * validation is required to infer "real" LOF. 
 * 
 * References: I used the LOF definition used in the 
 * following paper "A Systematic Survey of Loss-of-Function 
 * Variants in Human Protein-Coding Genes", Science, 2012
 * 
 * @author pcingola
 */
public class LossOfFunction {

	public LossOfFunction(Genome genome) {
	}

	/**
	 * Is this single change a LOF?
	 * @param changeEffect
	 * @return
	 */
	public boolean isLof(ChangeEffect changeEffect) {
		// The following effect types can be considered LOF
		switch (changeEffect.getEffectType()) {
		case SPLICE_SITE_ACCEPTOR:
		case SPLICE_SITE_DONOR:
			// Core splice sites are considered LOF
			if ((changeEffect.getMarker() != null) && (changeEffect.getMarker() instanceof SpliceSite)) {
				// Get splice site marker and check if it is 'core'
				SpliceSite spliceSite = (SpliceSite) changeEffect.getMarker();
				if (spliceSite.isCoreSpliceSite()) return true;

			}
			break;

		case STOP_GAINED:
			return true;

		case EXON_DELETED:
			return true;

		case FRAME_SHIFT:
			return true;

		case RARE_AMINO_ACID:
			// This one is not in the referenced papers, but we can assume that RARE AA changes are damaging.
			return true;

		case UTR_5_DELETED:
		case UTR_3_DELETED:
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
		}

		return lofCount > 0;
	}
}
