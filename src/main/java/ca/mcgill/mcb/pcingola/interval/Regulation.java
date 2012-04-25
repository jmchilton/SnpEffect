package ca.mcgill.mcb.pcingola.interval;

import java.util.List;

import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;

/**
 * Regulatory elements
 * 
 * @author pablocingolani
 */
public class Regulation extends Marker {

	private static final long serialVersionUID = -5607588295343642199L;

	String name = "";
	String cellType = "";

	public Regulation(Marker parent, int start, int end, int strand, String id, String name, String cellType) {
		super(parent, start, end, strand, id);
		type = EffectType.REGULATION.toString();
		this.name = name;
		this.cellType = cellType;
	}

	public String getCellType() {
		return cellType;
	}

	public String getName() {
		return name;
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
		EffectType effType = EffectType.REGULATION;
		changeEffect.set(this, effType, "");
		return changeEffect.newList();
	}

	@Override
	public String toString() {
		return getChromosomeName() + "\t" + start + "-" + end //
				+ " " //
				+ type + ((name != null) && (!name.isEmpty()) ? " '" + name + "'" : "");
	}

}
