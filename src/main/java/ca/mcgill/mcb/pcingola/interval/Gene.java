package ca.mcgill.mcb.pcingola.interval;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;
import ca.mcgill.mcb.pcingola.stats.ObservedOverExpectedCpG;

/**
 * Interval for a gene, as well as transcripts
 * 
 * @author pcingola
 *
 */
public class Gene extends IntervalAndSubIntervals<Transcript> implements Serializable {

	public enum GeneType {
		CODING, NON_CODING, UNKNOWN
	}

	private static final long serialVersionUID = 8419206759034068147L;

	String geneName;
	String bioType;

	public Gene(Marker parent, int start, int end, int strand, String id, String geneName, String bioType) {
		super(parent, start, end, strand, id);
		this.geneName = geneName;
		this.bioType = bioType;
		this.strand = strand;
		type = EffectType.GENE.toString();
	}

	/**
	 * Adjust start, end and strand values
	 * @return true if any adjustment was done
	 */
	public boolean adjust() {
		boolean changed = false;
		int strandSumGene = 0;
		int newStart = start, newEnd = start;

		if (newStart == 0 && newEnd == 0) {
			newStart = Integer.MAX_VALUE;
			newEnd = Integer.MIN_VALUE;
		}

		for (Transcript tr : this) {
			newStart = Math.min(newStart, tr.getStart());
			newEnd = Math.max(newEnd, tr.getEnd());

			for (Exon exon : tr.sortedStrand()) {
				newStart = Math.min(newStart, exon.getStart());
				newEnd = Math.max(newEnd, exon.getEnd());
				strandSumGene += exon.getStrand(); // Some exons have incorrect strands, we use the strand indicated by most exons
			}

			for (Utr utr : tr.getUtrs()) {
				newStart = Math.min(newStart, utr.getStart());
				newEnd = Math.max(newEnd, utr.getEnd());
			}
		}

		// Change gene strand?
		int newStrand = strandSumGene >= 0 ? 1 : -1;
		if (strand != newStrand) {
			strand = newStrand;
			changed = true;
		}

		// Change start?
		if (start != newStart) {
			start = newStart;
			changed = true;
		}

		// Change end?
		if (end != newEnd) {
			end = newEnd;
			changed = true;
		}

		return changed;
	}

	/**
	 * Calculate CpG bias: number of CpG / expected[CpG]
	 * @return
	 */
	public double cpgExonBias() {
		ObservedOverExpectedCpG oe = new ObservedOverExpectedCpG();
		return oe.oe(this);
	}

	public GeneType geneType() {
		if (bioType.length() > 0) {
			// Is it 'protein_coding'or a 'mRNA'?
			if (bioType.equalsIgnoreCase("protein_coding") || bioType.equalsIgnoreCase("mRNA")) return GeneType.CODING;
			return GeneType.NON_CODING;
		}

		return GeneType.UNKNOWN;
	}

	public String getBioType() {
		return bioType;
	}

	public String getGeneName() {
		return geneName;
	}

	@Override
	public int getStrand() {
		return strand;
	}

	/**
	 * Is any of the transcripts prtein coding?
	 * @return
	 */
	public boolean isProteinCoding() {
		for (Transcript tr : this)
			if (tr.isProteinCoding()) return true;
		return false;
	}

	@Override
	protected boolean isShowWarningIfParentDoesNotInclude() {
		return false;
	}

	/**
	 * Remove all non-canonical transcripts
	 */
	public void removeNonCanonical() {
		// Find canonical transcript (longest CDS)
		ArrayList<Transcript> toDelete = new ArrayList<Transcript>();
		Transcript canonical = null;
		for (Transcript t : this) {
			if ((canonical == null) || (canonical.cds().length() < t.cds().length())) canonical = t;
			toDelete.add(t);
		}

		// Found canonical? => Remove all others
		if (canonical != null) {
			// Remove all other transcripts
			toDelete.remove(canonical); // Do not remove canonical transcript.
			for (Transcript t : toDelete)
				remove(t);
		}
	}

	/**
	 * Get some details about the effect on this gene
	 * @param seqChange
	 * @return
	 */
	@Override
	public List<ChangeEffect> seqChangeEffect(SeqChange seqChange, ChangeEffect changeEffect) {
		if (!intersects(seqChange)) return ChangeEffect.emptyResults(); // Sanity check

		changeEffect.set(this, EffectType.GENE, "");

		boolean hitTranscript = false;

		ArrayList<ChangeEffect> changeEffectList = new ArrayList<ChangeEffect>();
		for (Transcript tr : this) {
			ChangeEffect chEff = changeEffect.clone();
			List<ChangeEffect> chEffList = tr.seqChangeEffect(seqChange, chEff);
			if (!chEffList.isEmpty()) {
				changeEffectList.addAll(chEffList);
				hitTranscript = true; // Did we hit any transcript?
			}
		}

		// May be none of the transcripts are actually hit 
		if (!hitTranscript) {
			changeEffect.set(this, EffectType.INTRAGENIC, "");
			return changeEffect.newList();
		}

		return changeEffectList;
	}

	public void setBioType(String bioType) {
		this.bioType = bioType;
	}

	/**
	 * Size of a genetic region for a given gene
	 * @param type
	 * @return
	 */
	public int sizeof(String type) {
		// Calculate size
		EffectType eff = EffectType.valueOf(type.toUpperCase());
		Markers all = new Markers();
		int len = 0;

		switch (eff) {
		case GENE:
			return size();

		case EXON:
			// Add all exons
			for (Transcript tr : this)
				for (Exon ex : tr)
					all.add(ex);
			break;

		case TRANSCRIPT:
			// Add all transcripts
			for (Transcript tr : this)
				all.add(tr);
			break;

		case INTRON:
			return Math.max(0, sizeof("TRANSCRIPT") - sizeof("EXON"));

		case UTR_3_PRIME:
			// Add all Utr3prime
			for (Transcript tr : this)
				for (Utr3prime utr : tr.get3primeUtrs())
					all.add(utr);
			break;

		case UTR_5_PRIME:
			// Add all Utr3prime
			for (Transcript tr : this)
				for (Utr5prime utr : tr.get5primeUtrs())
					all.add(utr);
			break;

		case UPSTREAM:
			for (Transcript tr : this)
				all.add(tr.getUpstream());

			break;

		case DOWNSTREAM:
			for (Transcript tr : this)
				all.add(tr.getDownstream());
			break;

		case SPLICE_SITE_ACCEPTOR:
			// Add all exons
			for (Transcript tr : this)
				for (Exon ex : tr)
					if (ex.getSpliceSiteAcceptor() != null) all.add(ex.getSpliceSiteAcceptor());
			break;

		case SPLICE_SITE_DONOR:
			// Add all exons
			for (Transcript tr : this)
				for (Exon ex : tr)
					if (ex.getSpliceSiteDonor() != null) all.add(ex.getSpliceSiteDonor());
			break;

		case INTRAGENIC:
			// We have to perform a set minus operation between this gene and all the transcripts
			Markers gene = new Markers();
			gene.add(this);

			// Create transcripts
			Markers trans = new Markers();
			for (Transcript tr : this)
				trans.add(tr);

			all = gene.minus(trans);
			break;

		case NONE:
			return 0;

		default:
			throw new RuntimeException("Unimplemented sizeof('" + type + "')");
		}

		// Merge and calculate total length
		Markers merged = all.merge();
		for (Marker m : merged)
			len += m.size();

		return len;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(getChromosomeName() + ":" + start + "-" + end);
		sb.append(", strand:" + strand);
		if ((id != null) && (id.length() > 0)) sb.append(", id:" + id);
		if ((geneName != null) && (geneName.length() > 0)) sb.append(", name:" + geneName);
		if ((bioType != null) && (bioType.length() > 0)) sb.append(", bioType:" + bioType);

		sb.append("\n");

		if (numChilds() > 0) {
			sb.append("Transcipts:\n");
			for (Transcript tint : sorted())
				sb.append("\t" + tint + "\n");
		}

		return sb.toString();
	}

}