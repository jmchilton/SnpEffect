package ca.mcgill.mcb.pcingola.interval;

import java.util.List;

import ca.mcgill.mcb.pcingola.motif.Pwm;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;

/**
 * Regulatory elements
 * 
 * @author pablocingolani
 */
public class Motif extends Marker {

	private static final long serialVersionUID = 8464487883781181867L;

	Pwm pwm;

	public Motif() {
		super();
		type = EffectType.MOTIF;
	}

	public Motif(Marker parent, int start, int end, int strand, String id, String name) {
		super(parent, start, end, strand, id);
		type = EffectType.MOTIF;
	}

	public Pwm getPwm() {
		return pwm;
	}

	/**
	 * Calculate the effect of this seqChange
	 * @param seqChange
	 * @param changeEffect
	 * @return
	 */
	@Override
	public List<ChangeEffect> seqChangeEffect(SeqChange seqChange, ChangeEffect changeEffect) {
		if (!intersects(seqChange)) return ChangeEffect.emptyResults(); // Sanity check
		EffectType effType = EffectType.MOTIF;
		changeEffect.set(this, effType, "");
		return changeEffect.newList();
	}

	public void setPwm(Pwm pwm) {
		this.pwm = pwm;
	}

}
