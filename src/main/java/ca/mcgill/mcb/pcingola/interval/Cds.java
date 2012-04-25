package ca.mcgill.mcb.pcingola.interval;

import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;

/**
 * Interval for a gene, as well as some other information: exons, utrs, cds, etc.
 * 
 * @author pcingola
 *
 */
public class Cds extends Marker {

	private static final long serialVersionUID = 1636197649250882952L;

	byte frame = 0;

	public Cds(Transcript parent, int start, int end, int strand, String id) {
		super(parent, start, end, strand, id);
		type = EffectType.CDS.toString();
	}

	public int getFrame() {
		return frame;
	}

	/**
	 * Frame can be {-1, 0, 1, 2}, where '-1' means unknown
	 * @param frame
	 */
	public void setFrame(int frame) {
		if( (frame > 2) || (frame < -1) ) throw new RuntimeException("Invalid frame value: " + frame);
		this.frame = (byte) frame;
	}

}
