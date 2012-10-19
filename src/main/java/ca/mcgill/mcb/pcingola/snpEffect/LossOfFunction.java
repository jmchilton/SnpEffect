package ca.mcgill.mcb.pcingola.snpEffect;

import java.util.List;

import ca.mcgill.mcb.pcingola.interval.Genome;

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
	 * Can this collection of effects produce a "Loss of function" 
	 * @param changeEffects
	 * @return
	 */
	public boolean isLof(List<ChangeEffect> changeEffects) {

		for (ChangeEffect changeEffect : changeEffects) {
			if (changeEffect.isSpliceSite()) {
				// TODO: We have to be carefull to define 'core' splice site (based on conservation?)

			}
			//			else if () {
			//
			//			}
		}

		return false;
	}

}
