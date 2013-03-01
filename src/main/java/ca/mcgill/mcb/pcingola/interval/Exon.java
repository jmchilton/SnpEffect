package ca.mcgill.mcb.pcingola.interval;

import ca.mcgill.mcb.pcingola.binseq.DnaNSequence;
import ca.mcgill.mcb.pcingola.binseq.DnaSequence;
import ca.mcgill.mcb.pcingola.interval.SeqChange.ChangeType;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.ErrorType;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.WarningType;
import ca.mcgill.mcb.pcingola.util.GprSeq;

/**
 * Interval for an exon
 * 
 * @author pcingola
 *
 */
public class Exon extends Marker {

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
	DnaSequence sequence;
	SpliceSiteAcceptor spliceSiteAcceptor;
	SpliceSiteDonor spliceSiteDonor;
	ExonSpliceType spliceType;

	protected Exon() {
		super();
		strand = 0;
		rank = 0;
		sequence = DnaSequence.empty();
		type = EffectType.EXON;
	}

	public Exon(Transcript parent, int start, int end, int strand, String id, int rank) {
		super(parent, start, end, strand, id);
		this.strand = (byte) strand;
		this.rank = rank;
		sequence = DnaSequence.empty();
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

		// Exon eliminated?
		if (ex == null) return null;

		// Sometimes 'apply' method return 'this'. Since we don't want to update the original exon, we have to create a clone
		if (ex == this) ex = (Exon) clone();

		// Update sites
		if (spliceSiteAcceptor != null) ex.spliceSiteAcceptor = (SpliceSiteAcceptor) spliceSiteAcceptor.apply(seqChange);
		if (spliceSiteDonor != null) ex.spliceSiteDonor = (SpliceSiteDonor) spliceSiteDonor.apply(seqChange);

		switch (seqChange.getChangeType()) {
		case SNP:
			applySnp(seqChange, ex);
			break;

		case INS:
			applyIns(seqChange, ex);
			break;

		case DEL:
			applyDel(seqChange, ex);
			break;

		case MNP:
			applyMnp(seqChange, ex);
			break;

		default:
			throw new RuntimeException("Unimplemented method for sequence change type " + seqChange.getChangeType());
		}

		return ex;
	}

	/**
	 * Apply a change type deletion
	 * @param seqChange
	 * @param ex
	 */
	protected void applyDel(SeqChange seqChange, Exon ex) {
		// Update sequence
		if ((sequence != null) && (!sequence.isEmpty())) {

			// Get sequence in positive strand direction
			String seq = isStrandPlus() ? sequence.getSequence() : sequence.reverseWc().getSequence();

			// Apply change to sequence
			int idxStart = seqChange.getStart() - start;
			int idxEnd = seqChange.getStart() - start + seqChange.size();

			StringBuilder newSeq = new StringBuilder();
			if (idxStart >= 0) newSeq.append(seq.substring(0, idxStart));
			if (idxEnd >= 0) newSeq.append(seq.substring(idxEnd));

			// Update sequence
			seq = newSeq.toString();
			ex.setSequence(isStrandPlus() ? seq : GprSeq.reverseWc(seq));
		}
	}

	/**
	 * Apply a change type insertion
	 * @param seqChange
	 * @param ex
	 */
	protected void applyIns(SeqChange seqChange, Exon ex) {
		// Update sequence
		if ((sequence != null) && (!sequence.isEmpty())) {

			// Get sequence in positive strand direction
			String seq = isStrandPlus() ? sequence.getSequence() : sequence.reverseWc().getSequence();

			String netChange = seqChange.netChange(this);
			// Apply change to sequence
			int idx = seqChange.getStart() - start - 1;
			if (idx >= 0) seq = seq.substring(0, idx + 1) + netChange + seq.substring(idx + 1);
			else seq = netChange + seq;

			// Update sequence
			ex.setSequence(isStrandPlus() ? seq : GprSeq.reverseWc(seq));
		}
	}

	/**
	 * Apply a change type MNP
	 * @param seqChange
	 * @param ex
	 */
	protected void applyMnp(SeqChange seqChange, Exon ex) {
		// Update sequence
		if ((sequence != null) && (!sequence.isEmpty())) {
			// Get sequence in positive strand direction
			String seq = isStrandPlus() ? sequence.getSequence() : sequence.reverseWc().getSequence();

			// Apply change to sequence
			int idxStart = seqChange.getStart() - start;
			int changeSize = seqChange.intersectSize(this);
			int idxEnd = seqChange.getStart() - start + changeSize;

			StringBuilder seqsb = new StringBuilder();
			seqsb.append(seq.substring(0, idxStart));
			seqsb.append(seqChange.getChange().substring(0, changeSize));
			seqsb.append(seq.substring(idxEnd));

			// Update sequence
			seq = seqsb.toString();
			ex.setSequence(isStrandPlus() ? seq : GprSeq.reverseWc(seq));
		}
	}

	/**
	 * Apply a change type SNP
	 * @param seqChange
	 * @param ex
	 */
	protected void applySnp(SeqChange seqChange, Exon ex) {
		// Update sequence
		if ((sequence != null) && (!sequence.isEmpty())) {
			// Get sequence in positive strand direction
			String seq = isStrandPlus() ? sequence.getSequence() : sequence.reverseWc().getSequence();

			// Apply change to sequence
			int idx = seqChange.getStart() - start;
			seq = seq.substring(0, idx) + seqChange.getChange() + seq.substring(idx + 1);

			// Update sequence
			ex.setSequence(isStrandPlus() ? seq : GprSeq.reverseWc(seq));
		}
	}

	/**
	 * Base in this exon at position 'index'
	 * @param idx
	 * @return
	 */
	public String basesAt(int index, int len) {
		if (strand < 0) {
			int idx = sequence.length() - index - len;
			return GprSeq.reverseWc(sequence.getBases(idx, len)); // Minus strand => Exon's sequence has been reversed and WC-complemented
		}
		return sequence.getBases(index, len);
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

	/**
	 * Get exon sequence
	 * 
	 * WARNING: Exon sequence is always according to coding 
	 * strand. E.g. if the strand is negative, the sequence 
	 * returned by this method is the reverse-WC that you see 
	 * in the reference genome
	 * 
	 * @param sequence
	 */
	public String getSequence() {
		return sequence.toString();
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

	/**
	 * Do we have a sequence for this exon?
	 * @return
	 */
	public boolean hasSequence() {
		if (size() <= 0) return true; // This interval has zero length, so sequence should be empty anyway (it is OK if its empty)
		return (sequence != null) && (!sequence.isEmpty());
	}

	@Override
	protected boolean isAdjustIfParentDoesNotInclude(Marker parent) {
		return true;
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

	/**
	 * Set exon sequence
	 * 
	 * WARNING: Exon sequence is always according to coding 
	 * strand. So use you should use setSequence( GprSeq.reverseWc( seq ) ) 
	 * if the exon is in negative strand.
	 * 
	 * @param sequence
	 */
	public void setSequence(String sequence) {
		if ((sequence == null) || (sequence.length() <= 0)) this.sequence = DnaSequence.empty();

		if (GprSeq.isAmbiguous(sequence)) this.sequence = new DnaNSequence(sequence); // Use DnaNSequence whch supports ambiguous sequences
		else this.sequence = new DnaSequence(sequence); // Use DnaSequence
	}

	@Override
	public String toString() {
		return getChromosomeName() + ":" + start + "-" + end //
				+ ((id != null) && (id.length() > 0) ? " '" + id + "'" : "") //
				+ " rank:" + rank //
				+ (sequence != null ? ", sequence: " + sequence : "");
	}

}
