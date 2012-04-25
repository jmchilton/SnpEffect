package ca.mcgill.mcb.pcingola.interval;

import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;

/**
 * Intron
 * 
 * @author pcingola
 *
 */
public class Intron extends Marker {

	private static final long serialVersionUID = -8283322526157264389L;

	public Intron(Transcript parent, int start, int end, int strand, String id) {
		super(parent, start, end, strand, id);
		type = EffectType.INTRON.toString();
	}

}
