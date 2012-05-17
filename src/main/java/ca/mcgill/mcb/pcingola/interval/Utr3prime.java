package ca.mcgill.mcb.pcingola.interval;

import java.util.List;

import ca.mcgill.mcb.pcingola.interval.SeqChange.ChangeType;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;

/**
 * Interval for a UTR (5 prime UTR and 3 prime UTR
 * 
 * @author pcingola
 *
 */
public class Utr3prime extends Utr {

	private static final long serialVersionUID = 5688641008301281991L;

	public Utr3prime(Marker parent, int start, int end, int strand, String id) {
		super(parent, start, end, strand, id);
		type = EffectType.UTR_3_PRIME.toString();
	}

	@Override
	public List<ChangeEffect> seqChangeEffect(SeqChange seqChange, ChangeEffect changeEffect) {
		if (seqChange.includes(this) && (seqChange.getChangeType() == ChangeType.DEL)) {
			changeEffect.set(this, EffectType.UTR_3_DELETED, ""); // A UTR was removed entirely
			return changeEffect.newList();
		}

		Transcript tint = (Transcript) findParent(Transcript.class);
		Exon exon = (Exon) findParent(Exon.class);

		String utrDistStr = utrDistance(seqChange, tint);

		changeEffect.set(this, EffectType.valueOf(type), utrDistStr);
		changeEffect.setExon(exon);

		Exon eint = (Exon) findParent(Exon.class);
		if (eint != null) eint.check(seqChange, changeEffect); // Check that base matches the expected one

		return changeEffect.newList();
	}

	/**
	 * Calculate distance from beginning of 3'UTRs
	 * 
	 * @param snp
	 * @param utr
	 * @return
	 */
	@Override
	String utrDistance(SeqChange snp, Transcript tint) {
		List<Utr3prime> utrs = tint.get3primeUtrs();
		boolean fromEnd = strand < 0; // We want distance from end of transcript (beginning of 3'UTR)
		int dist = snp.distanceFrom(utrs, fromEnd) + 1;
		return dist + " bases from CDS";
	}

}