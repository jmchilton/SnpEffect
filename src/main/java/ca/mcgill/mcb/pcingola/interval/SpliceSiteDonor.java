package ca.mcgill.mcb.pcingola.interval;

import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;

/**
 * Interval for a splice site donnor
 * 
 * Note: Splice sites donnor are defined as the first 2 bases of an intron 
 * Reference: http://en.wikipedia.org/wiki/RNA_splicing
 * 
 * @author pcingola
 *
 */
public class SpliceSiteDonor extends SpliceSite {

	private static final long serialVersionUID = -2117470153797320999L;

	public SpliceSiteDonor() {
		super();
		type = EffectType.SPLICE_SITE_DONOR;
	}

	public SpliceSiteDonor(Marker parent, int start, int end, int strand, String id) {
		super(parent, start, end, strand, id);
		type = EffectType.SPLICE_SITE_DONOR;
	}
}
