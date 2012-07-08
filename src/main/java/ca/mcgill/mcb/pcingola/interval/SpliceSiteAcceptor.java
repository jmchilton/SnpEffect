package ca.mcgill.mcb.pcingola.interval;

import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;

/**
 * Interval for a splice site acceptor
 * 
 * Note: Splice sites donnor are defined as the last 2 bases of an intron 
 * Reference: http://en.wikipedia.org/wiki/RNA_splicing
 * 
 * @author pcingola
 *
 */
public class SpliceSiteAcceptor extends SpliceSite {

	private static final long serialVersionUID = -7416687954435361328L;

	public SpliceSiteAcceptor() {
		super();
		type = EffectType.SPLICE_SITE_ACCEPTOR;
	}

	public SpliceSiteAcceptor(Marker parent, int start, int end, int strand, String id) {
		super(parent, start, end, strand, id);
		type = EffectType.SPLICE_SITE_ACCEPTOR;
	}
}
