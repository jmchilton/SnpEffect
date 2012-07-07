package ca.mcgill.mcb.pcingola.interval;

import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;

public class RareAminoAcid extends Marker {

	private static final long serialVersionUID = -1926572865764543849L;

	public RareAminoAcid(Marker parent, int start, int end, String id) {
		super(parent, start, end, 1, id);
		type = EffectType.RARE_AMINO_ACID;
	}
}
