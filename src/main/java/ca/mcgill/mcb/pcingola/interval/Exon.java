package ca.mcgill.mcb.pcingola.interval;

import ca.mcgill.mcb.pcingola.interval.SeqChange.ChangeType;
import ca.mcgill.mcb.pcingola.serializer.MarkerSerializer;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.ErrorType;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.WarningType;

/**
 * Interval for an exon
 * 
 * @author pcingola
 *
 */
public class Exon extends MarkerSeq {

	/**
	 * Characterize exons based on alternative splicing
	 * References: "Alternative splicing and evolution - diversification, exon definition and function"  (see Box 1)
	 */
	public enum ExonSpliceType {
		RETAINED, // All transcripts have this exon
		SKIPPED, // Some transcripts skip it
		ALTTENATIVE_3SS, // Some transcripts have and alternative 3' exon start 
		ALTTENATIVE_5SS, // Some transcripts have and alternative 5' exon end
		MUTUALLY_EXCLUSIVE, // Mutually exclusive (respect to other exon)
		ALTTENATIVE_PROMOMOTER, // The first exon is different in some transcripts.
		ALTTENATIVE_POLY_A, // The last exon.
	}

	private static final long serialVersionUID = 5324352193278472543L;

	byte frame = 0;
	int rank; // Exon rank in transcript
	SpliceSiteAcceptor spliceSiteAcceptor;
	SpliceSiteDonor spliceSiteDonor;
	ExonSpliceType spliceType;

	public Exon() {
		super();
		rank = 0;
		type = EffectType.EXON;
	}

	public Exon(Transcript parent, int start, int end, int strand, String id, int rank) {
		super(parent, start, end, strand, id);
		this.strand = (byte) strand;
		this.rank = rank;
		type = EffectType.EXON;
	}

	/**
	 * Apply seqChange to exon
	 * 
	 * WARNING: There might be conditions which change the exon type (e.g. an intron is deleted)
	 * 			Nevertheless ExonSpliceType s not updated since it reflects the exon type before a sequence change. 
	 * 
	 */
	@Override
	public Exon apply(SeqChange seqChange) {
		// Create new exon with updated coordinates
		Exon ex = (Exon) super.apply(seqChange);

		// Update sites
		if (spliceSiteAcceptor != null) ex.spliceSiteAcceptor = (SpliceSiteAcceptor) spliceSiteAcceptor.apply(seqChange);
		if (spliceSiteDonor != null) ex.spliceSiteDonor = (SpliceSiteDonor) spliceSiteDonor.apply(seqChange);

		return ex;
	}

	/**
	 * Check that the base in the exon corresponds with the one in the SNP
	 * @param seqChange
	 * @param results
	 */
	public void check(SeqChange seqChange, ChangeEffect results) {
		// Only makes sense for SNPs and MNPs
		if ((seqChange.getChangeType() != ChangeType.SNP) && (seqChange.getChangeType() != ChangeType.MNP)) return;

		int mstart = Math.max(seqChange.getStart(), start);
		int idxStart = mstart - start;
		if (sequence.length() <= 0) results.addWarning(WarningType.WARNING_SEQUENCE_NOT_AVAILABLE);
		else if (idxStart >= sequence.length()) results.addError(ErrorType.ERROR_OUT_OF_EXON);
		else {
			int mend = Math.min(seqChange.getEnd(), end);
			int len = mend - mstart + 1;

			String realReference = basesAt(idxStart, len).toUpperCase();
			String changeReference = seqChange.reference().substring(mstart - seqChange.getStart(), mend - seqChange.getStart() + 1);

			// Reference sequence different than expected?
			if (!realReference.equals(changeReference)) results.addWarning(WarningType.WARNING_REF_DOES_NOT_MATCH_GENOME);
		}
	}

	/**
	 * Create a splice site acceptor of 'maxSize' length
	 * @param size
	 * @return
	 */
	public SpliceSiteAcceptor createSpliceSiteAcceptor(int size) {
		size = size - 1;
		if (size < 0) return null;

		if (strand >= 0) spliceSiteAcceptor = new SpliceSiteAcceptor(this, start - 1 - size, start - 1, strand, id);
		else spliceSiteAcceptor = new SpliceSiteAcceptor(this, end + 1, end + 1 + size, strand, id);

		return spliceSiteAcceptor;
	}

	/**
	 * Create a splice site donor of 'maxSize' length
	 * @param size
	 * @return
	 */
	public SpliceSiteDonor createSpliceSiteDonor(int size) {
		size = size - 1;
		if (size < 0) return null;

		if (strand >= 0) spliceSiteDonor = new SpliceSiteDonor(this, end + 1, end + 1 + size, strand, id);
		else spliceSiteDonor = new SpliceSiteDonor(this, start - 1 - size, start - 1, strand, id);

		return spliceSiteDonor;
	}

	public int getFrame() {
		return frame;
	}

	public int getRank() {
		return rank;
	}

	public SpliceSiteAcceptor getSpliceSiteAcceptor() {
		return spliceSiteAcceptor;
	}

	public SpliceSiteDonor getSpliceSiteDonor() {
		return spliceSiteDonor;
	}

	public ExonSpliceType getSpliceType() {
		return spliceType;
	}

	@Override
	protected boolean isAdjustIfParentDoesNotInclude(Marker parent) {
		return true;
	}

	/**
	 * Query all genomic regions that intersect 'marker'
	 */
	public Markers query(Marker marker) {
		Markers markers = new Markers();
		if (marker.intersects(spliceSiteAcceptor)) markers.add(spliceSiteAcceptor);
		if (marker.intersects(spliceSiteDonor)) markers.add(spliceSiteDonor);
		return markers;
	}

	/**
	 * Parse a line from a serialized file
	 * @param line
	 * @return
	 */
	@Override
	public void serializeParse(MarkerSerializer markerSerializer) {
		super.serializeParse(markerSerializer);
		frame = (byte) markerSerializer.getNextFieldInt();
		rank = markerSerializer.getNextFieldInt();
		setSequence(markerSerializer.getNextField());
		spliceSiteDonor = (SpliceSiteDonor) markerSerializer.getNextFieldMarker();
		spliceSiteAcceptor = (SpliceSiteAcceptor) markerSerializer.getNextFieldMarker();

		String exType = markerSerializer.getNextField();
		if ((exType != null) && !exType.isEmpty()) spliceType = ExonSpliceType.valueOf(exType);
	}

	/**
	 * Create a string to serialize to a file
	 * @return
	 */
	@Override
	public String serializeSave(MarkerSerializer markerSerializer) {
		int ssdId = markerSerializer.save(spliceSiteDonor);
		int ssaId = markerSerializer.save(spliceSiteAcceptor);

		return super.serializeSave(markerSerializer) //
				+ "\t" + frame //
				+ "\t" + rank //
				+ "\t" + sequence //
				+ "\t" + ssdId //
				+ "\t" + ssaId //
				+ "\t" + (spliceType != null ? spliceType.toString() : "")//				
		;
	}

	/**
	 * Frame can be {-1, 0, 1, 2}, where '-1' means unknown
	 * @param frame
	 */
	public void setFrame(int frame) {
		if ((frame > 2) || (frame < -1)) throw new RuntimeException("Invalid frame value: " + frame);
		this.frame = (byte) frame;
	}

	public void setRank(int rank) {
		this.rank = rank;
	}

	@Override
	public String toString() {
		return getChromosomeName() + ":" + start + "-" + end //
				+ ((id != null) && (id.length() > 0) ? " '" + id + "'" : "") //
				+ " rank:" + rank //
				+ (sequence != null ? ", sequence: " + sequence : "");
	}

}
