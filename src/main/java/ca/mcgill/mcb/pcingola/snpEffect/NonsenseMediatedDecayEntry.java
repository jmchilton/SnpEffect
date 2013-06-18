package ca.mcgill.mcb.pcingola.snpEffect;

import ca.mcgill.mcb.pcingola.interval.Gene;

/**
 * A Nonsense Mediated Decay (NMD) entry in a VCF file  
 * 
 * @author pcingola
 */
public class NonsenseMediatedDecayEntry extends LossOfFunctionEntry {

	public NonsenseMediatedDecayEntry(Gene gene, double percentOfTranscriptsAffected) {
		super(gene, percentOfTranscriptsAffected);
	}

	public NonsenseMediatedDecayEntry(String vcfNmdStr) {
		super(vcfNmdStr);
	}

	public NonsenseMediatedDecayEntry(String geneName, String geneId, int transcriptInGene, double percentOfTranscriptsAffected) {
		super(geneName, geneId, transcriptInGene, percentOfTranscriptsAffected);
	}

}
