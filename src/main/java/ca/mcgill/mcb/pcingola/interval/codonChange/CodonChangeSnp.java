package ca.mcgill.mcb.pcingola.interval.codonChange;

import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;

/**
 * Calculate codon changes produced by a SNP
 * @author pcingola
 */
public class CodonChangeSnp extends CodonChange {

	public CodonChangeSnp(SeqChange seqChange, Transcript transcript, ChangeEffect changeEffect) {
		super(seqChange, transcript, changeEffect);
		returnNow = true; // A SNP can only affect one exon
	}

	/**
	 * Analyze SNPs in this transcript.
	 * Add changeEffect to 'changeEffect'
	 */
	@Override
	boolean codonChangeSingle(ChangeEffect changeEffect, Exon exon) {
		// Get old and new codons
		codonsOld = codonsOld();
		codonsNew = codonsNew();
		changeEffect.set(transcript, EffectType.CODON_CHANGE, "");
		changeEffect.setCodons(codonsOld, codonsNew, codonNum, codonIndex);
		return true;
	}

	/**
	 * Get new (modified) codons 
	 * @return
	 */
	@Override
	String codonsNew() {
		char codonChars[] = codonsOld.toLowerCase().toCharArray();
		char snpBase = seqChange.netChange(transcript.getStrand()).charAt(0);
		codonChars[codonIndex] = Character.toUpperCase(snpBase);

		String codonsNew = new String(codonChars);
		return codonsNew;
	}

	/**
	 * Get original codons in CDS
	 * @param codonNum
	 * @return
	 */
	@Override
	public String codonsOld() {
		int numCodons = 1;
		int minBase = codonNum * CodonChange.CODON_SIZE;
		if (minBase < 0) minBase = 0;

		int len = transcript.cds().length();
		int maxBase = codonNum * CodonChange.CODON_SIZE + numCodons * CodonChange.CODON_SIZE;
		if (maxBase > len) maxBase = len;

		char codonChars[] = transcript.cds().substring(minBase, maxBase).toLowerCase().toCharArray();
		codonChars[codonIndex] = Character.toUpperCase(codonChars[codonIndex]);
		String codon = new String(codonChars);
		return codon;
	}
}
