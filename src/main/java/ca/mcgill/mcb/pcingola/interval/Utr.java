package ca.mcgill.mcb.pcingola.interval;

/**
 * Interval for a UTR (5 prime UTR and 3 prime UTR
 * 
 * @author pcingola
 *
 */
public abstract class Utr extends Marker {

	private static final long serialVersionUID = 1636197649250882952L;

	public Utr() {
		super();
	}

	public Utr(Marker parent, int start, int end, int strand, String id) {
		super(parent, start, end, strand, id);
	}

	abstract String utrDistance(SeqChange snp, Transcript tint);

}
