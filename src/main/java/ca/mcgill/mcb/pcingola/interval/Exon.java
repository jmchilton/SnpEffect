package ca.mcgill.mcb.pcingola.interval;

import java.util.List;

import ca.mcgill.mcb.pcingola.binseq.DnaNSequence;
import ca.mcgill.mcb.pcingola.binseq.DnaSequence;
import ca.mcgill.mcb.pcingola.interval.SeqChange.ChangeType;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;
import ca.mcgill.mcb.pcingola.util.GprSeq;

/**
 * Interval for an exon
 * 
 * @author pcingola
 *
 */
public class Exon extends Marker {

	/**
	 * Caracterize exons based on alternative splicing
	 * References: "Alternative splicing and evolution - diversification, exon definition and function"  (see Box 1)
	 */
	public enum ExonSpliceType {
		RETAINED, // All transcripts have it
		SKIPPED, // Some transcripts skip it
		ALTTENATIVE_3SS, // Some transcripts have and alternative 3' exon start 
		ALTTENATIVE_5SS, // Some transcripts have and alternative 5' exon end
		MUTUALLY_EXCLUSIVE, // Mutually exclusive (respect to other exon)
		ALTTENATIVE_PROMOMOTER, // The first exon is different in some transcripts.
		ALTTENATIVE_POLY_A, // The last exon.
	}

	private static final long serialVersionUID = 5324352193278472543L;
	public static final int SPLICE_SITE_SIZE = 2;;

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

	public Exon(Marker parent, int start, int end, int strand, String id, int rank) {
		super(parent, start, end, strand, id);
		this.strand = (byte) strand;
		this.rank = rank;
		sequence = DnaSequence.empty();
		type = EffectType.EXON;
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
		if (sequence.length() <= 0) {
			results.addWarning("WARNING_SEQUENCE_NOT_AVAILABLE");
		} else if (idxStart >= sequence.length()) {
			results.addError("ERROR_OUT_OF_EXON");
		} else {
			int mend = Math.min(seqChange.getEnd(), end);
			int len = mend - mstart + 1;

			String realReference = basesAt(idxStart, len).toUpperCase();
			String changeReference = seqChange.reference().substring(mstart - seqChange.getStart(), mend - seqChange.getStart() + 1);

			// Reference sequence different than expected?
			if (!realReference.equals(changeReference)) {
				results.addWarning("WARNING_REF_IS_" + realReference + "_NOT_" + seqChange.getReference());
			}
		}
	}

	public SpliceSiteAcceptor createSpliceSiteAcceptor(int maxSize) {
		maxSize = Math.min(SPLICE_SITE_SIZE, maxSize) - 1;
		if (maxSize < 0) return null;

		if (strand >= 0) spliceSiteAcceptor = new SpliceSiteAcceptor(this, start - 1 - maxSize, start - 1, strand, id);
		else spliceSiteAcceptor = new SpliceSiteAcceptor(this, end + 1, end + 1 + maxSize, strand, id);

		return spliceSiteAcceptor;
	}

	public SpliceSiteDonor createSpliceSiteDonor(int maxSize) {
		maxSize = Math.min(SPLICE_SITE_SIZE, maxSize) - 1;
		if (maxSize < 0) return null;

		if (strand >= 0) spliceSiteDonor = new SpliceSiteDonor(this, end + 1, end + 1 + maxSize, strand, id);
		else spliceSiteDonor = new SpliceSiteDonor(this, start - 1 - maxSize, start - 1, strand, id);

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

	@Override
	protected boolean isAdjustIfParentDoesNotInclude(Marker parent) {
		return true;
	}

	/**
	 * Get some details about the effect on this transcript
	 * @param seqChange
	 * @return
	 */
	@Override
	public List<ChangeEffect> seqChangeEffect(SeqChange seqChange, ChangeEffect changeEffect) {
		if (!intersects(seqChange)) return ChangeEffect.emptyResults(); // Sanity check

		// An exon was removed entirely?
		if (seqChange.includes(this) && seqChange.isDel()) {
			changeEffect.setCodons("", "", -1, -1);
			changeEffect.set(this, EffectType.EXON_DELETED, "");
			return changeEffect.newList();
		}

		changeEffect.setExon(this); // We do NOT use changeEffect.set(this, EXON, "") because the information has already been set by 'Transcript'
		check(seqChange, changeEffect); // Check that the base in the exon corresponds with the one in the SNP
		return changeEffect.newList();
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
