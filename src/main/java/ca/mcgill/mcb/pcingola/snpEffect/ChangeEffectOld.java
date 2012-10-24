package ca.mcgill.mcb.pcingola.snpEffect;

import java.util.ArrayList;
import java.util.List;

import ca.mcgill.mcb.pcingola.codons.CodonTable;
import ca.mcgill.mcb.pcingola.interval.Custom;
import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Regulation;
import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.interval.Transcript;

/**
 * Object used to show results of a sequence change effect.
 * 
 * TODO: Several things to improve in this class 
 * 	- There should be a class called 'ChangeEffects' that accumulated multiple effects from a seqChange
 * 	- 
 * 
 * 
 * @author pcingola
 */
public class ChangeEffectOld implements Cloneable {

	public enum Coding {
		CODING, NON_CODING
	}

	public enum EffectImpact {
		HIGH, MODERATE, LOW, MODIFIER
	}

	public enum EffectType {
		NONE //
		, CHROMOSOME //
		, INTERGENIC //
		, UPSTREAM //
		, UTR_5_PRIME //
		, UTR_5_DELETED //		
		, START_GAINED //
		, SPLICE_SITE_ACCEPTOR //
		, SPLICE_SITE_BRANCH //
		, SPLICE_SITE_BRANCH_U12 //
		, SPLICE_SITE_DONOR //
		, START_LOST // 
		, SYNONYMOUS_START //
		, NON_SYNONYMOUS_START //
		, CDS //
		, GENE //
		, GENOME //
		, TRANSCRIPT //
		, EXON //		
		, EXON_DELETED //		
		, NON_SYNONYMOUS_CODING //
		, SYNONYMOUS_CODING //
		, FRAME_SHIFT //
		, CODON_CHANGE //
		, CODON_INSERTION //
		, CODON_CHANGE_PLUS_CODON_INSERTION //
		, CODON_DELETION //
		, CODON_CHANGE_PLUS_CODON_DELETION //
		, RARE_AMINO_ACID //
		, STOP_GAINED //
		, SYNONYMOUS_STOP //
		, NON_SYNONYMOUS_STOP //
		, STOP_LOST //
		, INTRON //
		, UTR_3_PRIME //
		, UTR_3_DELETED //		
		, DOWNSTREAM //
		, INTRON_CONSERVED //
		, INTERGENIC_CONSERVED //
		, INTRAGENIC //
		, REGULATION //
		, MICRO_RNA //
		, CUSTOM
	};

	/**
	 * This class is only getFused for SNPs
	 */
	public enum FunctionalClass {
		NONE, SILENT, MISSENSE, NONSENSE
	};

	static final boolean COMPATIBLE_v1_8 = true; // Activate this in order to get the same out as version 1.8. This is only for testing & debugging 

	SeqChange seqChange = null;
	EffectType effectType = EffectType.NONE;
	EffectImpact effectImpact = null;
	Marker marker = null;
	Exon exon = null;
	String error = "", warning = "", message = ""; // Any message, warning or error?
	String codonsOld = "", codonsNew = ""; // Codon change information
	String codonsAroundOld = "", codonsAroundNew = ""; // Codons arround
	int codonNum = -1; // Codon number
	int codonIndex = -1; // Index within a codon
	int codonDegeneracy = -1; // Codon degeneracy
	String aaOld = "", aaNew = ""; // Amino acid changes
	String aasAroundOld = "", aasAroundNew = ""; // Amino acids arround

	/**
	 *  An empty list of results;
	 * @return
	 */
	public static List<ChangeEffectOld> emptyResults() {
		return new ArrayList<ChangeEffectOld>();
	}

	public ChangeEffectOld(SeqChange seqChange) {
		this.seqChange = seqChange;
	}

	public void addError(String err) {
		error += (error.isEmpty() ? "" : "|") + err;
	}

	public void addWarning(String warn) {
		warning += (warning.isEmpty() ? "" : "|") + warn;
	}

	@Override
	public ChangeEffectOld clone() {
		try {
			return (ChangeEffectOld) super.clone();
		} catch (CloneNotSupportedException e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * Create a string for codon effect
	 * @param showAaChange : If true, include codon change, biotype, etc.
	 * @return
	 */
	String codonEffect(boolean showAaChange, boolean showBioType) {
		String codonEffect = "";
		if ((marker == null) || (codonNum < 0)) return codonEffect;

		// Add codon effect
		codonEffect += effectType;

		// Append codon change
		if (showAaChange) codonEffect += "(" + getAaChange() + ")";

		return codonEffect;
	}

	/**
	 * Show a string with overall effect
	 * @param shortFormat
	 * @param showAaChange
	 * @return
	 */
	public String effect(boolean shortFormat, boolean showAaChange, boolean showBioType) {
		String e = "";
		String codonEffect = codonEffect(showAaChange, showBioType); // Codon effect

		if (!codonEffect.isEmpty()) e = codonEffect;
		else if (isRegulation()) e = effectType.toString() + "[" + ((Regulation) marker).getName() + "]";
		else if (isIntergenic() || isIntron() || isSpliceSite()) e = effectType.toString();
		else if (!message.isEmpty()) e = effectType.toString() + ": " + message;
		else if (marker == null) e = effectType.toString(); // There are cases when no marker is associated (e.g. "Out of chromosome", "No such chromosome", etc.)
		else e = effectType.toString() + ": " + marker.getId();

		if (shortFormat) e = e.split(":")[0];
		return e;
	}

	/**
	 * Amino acid change string
	 * @return
	 */
	public String getAaChange() {
		if (aaOld.isEmpty() && aaNew.isEmpty()) return "";
		if (aaOld.equals(aaNew)) return aaNew;
		return aaOld + "/" + aaNew;
	}

	/**
	 * Amino acid change string (HGVS style)
	 * References: http://www.hgvs.org/mutnomen/recs.html
	 * 
	 * @return
	 */
	public String getAaChangeHgvs() {
		if (aaOld.isEmpty() && aaNew.isEmpty()) {
			if (codonNum >= 0) return "" + (codonNum + 1);
			return "";
		}
		if (aaOld.equals(aaNew)) return aaNew + (codonNum + 1);
		return aaOld + (codonNum + 1) + aaNew;
	}

	/**
	 * Amino acid length (negative if there is none)
	 * @return Amino acid length (CDS length / 3 ) or '-1' if there is no CDS length
	 */
	public int getAaLength() {
		int cdsLen = getCdsLength();
		if (cdsLen < 0) return -1;

		int lenNoStop = Math.max(0, cdsLen - 3); // Do not include the STOP codon
		return lenNoStop / 3;
	}

	public String getAaNew() {
		return aaNew;
	}

	public String getAaOld() {
		return aaOld;
	}

	/**
	 * CDS length (negative if there is none)
	 * @return
	 */
	public int getCdsLength() {
		if (marker == null) return -1;

		// CDS size info
		Transcript tr = null;
		if (exon != null) tr = (Transcript) exon.findParent(Transcript.class);
		else tr = (Transcript) marker.findParent(Transcript.class);

		if ((tr != null) && tr.isProteinCoding()) return tr.cds().length();
		return -1;
		//		int cdsSize = -1;
		//
		//		if (COMPATIBLE_v1_8 && (isUpstream() || isDownstream() || isUtr() || isSpliceSite() || isStartGained())) cdsSize = -1;
		//		return cdsSize;
	}

	/**
	 * Codon change string
	 * @return
	 */
	public String getCodonChange() {
		if (codonsOld.isEmpty() && codonsNew.isEmpty()) return "";
		if (codonsOld.equals(codonsNew)) return codonsNew;
		return codonsOld + "/" + codonsNew;
	}

	public String getCodonsNew() {
		return codonsNew;
	}

	public String getCodonsOld() {
		return codonsOld;
	}

	/**
	 * Return impact of this effect 
	 * @return
	 */
	public EffectImpact getEffectImpact() {
		if (effectImpact == null) {
			switch (effectType) {
			case EXON_DELETED:
			case FRAME_SHIFT:
			case SPLICE_SITE_ACCEPTOR:
			case SPLICE_SITE_DONOR:
			case START_LOST:
			case STOP_GAINED:
			case STOP_LOST:
			case RARE_AMINO_ACID:
				effectImpact = EffectImpact.HIGH;
				break;

			case CODON_CHANGE:
			case CODON_CHANGE_PLUS_CODON_DELETION:
			case CODON_CHANGE_PLUS_CODON_INSERTION:
			case CODON_DELETION:
			case CODON_INSERTION:
			case NON_SYNONYMOUS_CODING:
			case SPLICE_SITE_BRANCH_U12:
			case UTR_3_DELETED:
			case UTR_5_DELETED:
				effectImpact = EffectImpact.MODERATE;
				break;

			case SPLICE_SITE_BRANCH:
			case NON_SYNONYMOUS_START:
			case NON_SYNONYMOUS_STOP:
			case START_GAINED:
			case SYNONYMOUS_CODING:
			case SYNONYMOUS_START:
			case SYNONYMOUS_STOP:
				effectImpact = EffectImpact.LOW;
				break;

			case CDS:
			case CHROMOSOME:
			case CUSTOM:
			case DOWNSTREAM:
			case EXON:
			case GENE:
			case INTRAGENIC:
			case INTERGENIC:
			case INTERGENIC_CONSERVED:
			case INTRON:
			case INTRON_CONSERVED:
			case NONE:
			case REGULATION:
			case TRANSCRIPT:
			case UPSTREAM:
			case UTR_3_PRIME:
			case UTR_5_PRIME:
				effectImpact = EffectImpact.MODIFIER;
				break;

			default:
				throw new RuntimeException("Unknown impact for effect type: '" + effectType + "'");
			}
		}
		return effectImpact;
	}

	public EffectType getEffectType() {
		return effectType;
	}

	public String getError() {
		return error;
	}

	public Exon getExon() {
		if (exon != null) return exon;

		// Try to find exon
		if (marker != null) {
			Exon e = (Exon) marker.findParent(Exon.class);
			if (e != null) return e;
		}

		return null;
	}

	/**
	 * Return functional class of this effect (i.e.  NONSENSE, MISSENSE, SILENT or NONE)
	 * @return
	 */
	public FunctionalClass getFunctionalClass() {
		if (seqChange.isSnp()) {
			if (!aaNew.equals(aaOld)) {
				CodonTable codonTable = marker.codonTable();
				if (codonTable.isStop(codonsNew)) return FunctionalClass.NONSENSE;

				return FunctionalClass.MISSENSE;
			}
			if (!codonsNew.equals(codonsOld)) return FunctionalClass.SILENT;
		}

		return FunctionalClass.NONE;
	}

	public Gene getGene() {
		if (marker == null) return null;
		return (Gene) marker.findParent(Gene.class);
	}

	public String getGeneRegion() {
		switch (effectType) {
		case NONE:
		case CHROMOSOME:
		case CUSTOM:
		case CDS:
			return EffectType.NONE.toString();
		case INTERGENIC:
		case INTERGENIC_CONSERVED:
			return EffectType.INTERGENIC.toString();
		case UPSTREAM:
			return EffectType.UPSTREAM.toString();
		case UTR_5_PRIME:
		case UTR_5_DELETED:
		case START_GAINED:
			return EffectType.UTR_5_PRIME.toString();
		case SPLICE_SITE_ACCEPTOR:
			return EffectType.SPLICE_SITE_ACCEPTOR.toString();
		case SPLICE_SITE_BRANCH_U12:
		case SPLICE_SITE_BRANCH:
			return EffectType.SPLICE_SITE_BRANCH.toString();
		case SPLICE_SITE_DONOR:
			return EffectType.SPLICE_SITE_DONOR.toString();
		case INTRAGENIC:
		case START_LOST:
		case SYNONYMOUS_START:
		case NON_SYNONYMOUS_START:
		case GENE:
		case TRANSCRIPT:
			if (hasExon()) return EffectType.EXON.toString();
			return EffectType.NONE.toString();
		case EXON:
		case EXON_DELETED:
		case NON_SYNONYMOUS_CODING:
		case SYNONYMOUS_CODING:
		case FRAME_SHIFT:
		case CODON_CHANGE:
		case CODON_INSERTION:
		case CODON_CHANGE_PLUS_CODON_INSERTION:
		case CODON_DELETION:
		case CODON_CHANGE_PLUS_CODON_DELETION:
		case STOP_GAINED:
		case SYNONYMOUS_STOP:
		case NON_SYNONYMOUS_STOP:
		case STOP_LOST:
		case RARE_AMINO_ACID:
			return EffectType.EXON.toString();
		case INTRON:
		case INTRON_CONSERVED:
			return EffectType.INTRON.toString();
		case UTR_3_PRIME:
		case UTR_3_DELETED:
			return EffectType.UTR_3_PRIME.toString();
		case DOWNSTREAM:
			return EffectType.DOWNSTREAM.toString();
		case REGULATION:
			return EffectType.REGULATION.toString();
		default:
			throw new RuntimeException("Unknown gene region for effect type: '" + effectType + "'");
		}
	}

	public Marker getMarker() {
		return marker;
	}

	public SeqChange getSeqChange() {
		return seqChange;
	}

	public Transcript getTranscript() {
		if (marker == null) return null;
		return (Transcript) marker.findParent(Transcript.class);
	}

	public String getWarning() {
		return warning;
	}

	public boolean hasExon() {
		return (exon != null);
	}

	public boolean hasWarning() {
		return (warning != null) && (!warning.isEmpty());
	}

	/**
	 * Show data header
	 * @return
	 */
	public String header() {
		return "Warnings\t" //
				+ "Gene_ID\t" //
				+ "Gene_name\t" //
				+ "Bio_type\t" //
				+ "Trancript_ID\t" //
				+ "Exon_ID\t" //
				+ "Exon_Rank\t" //
				+ "Effect\t" //
				+ "old_AA/new_AA\t" //
				+ "Old_codon/New_codon\t" //
				+ "Codon_Num(CDS)\t" //
				+ "Codon_Degeneracy\t" //
				+ "CDS_size\t" //
				+ "Codons_around\t" //
				+ "AAs_around\t" //
				+ "Custom_interval_ID";
	}

	public boolean isCustom() {
		return effectType == EffectType.CUSTOM;
	}

	public boolean isDownstream() {
		return effectType == EffectType.DOWNSTREAM;
	}

	public boolean isExon() {
		return hasExon() || (effectType == EffectType.EXON_DELETED);
	}

	public boolean isFrameShift() {
		return (effectType == EffectType.FRAME_SHIFT);
	}

	public boolean isIntergenic() {
		return (effectType == EffectType.INTERGENIC) || (effectType == EffectType.INTERGENIC_CONSERVED);
	}

	public boolean isIntron() {
		return effectType == EffectType.INTRON;
	}

	public boolean isRegulation() {
		return (effectType == EffectType.REGULATION);
	}

	public boolean isSpliceSite() {
		return (effectType == EffectType.SPLICE_SITE_DONOR) //
				|| (effectType == EffectType.SPLICE_SITE_ACCEPTOR) //
				|| (effectType == EffectType.SPLICE_SITE_BRANCH) //
				|| (effectType == EffectType.SPLICE_SITE_BRANCH_U12) //
		;
	}

	public boolean isStartGained() {
		return effectType == EffectType.START_GAINED;
	}

	public boolean isUpstream() {
		return (effectType == EffectType.UPSTREAM) || (effectType == EffectType.START_GAINED);
	}

	public boolean isUtr() {
		return (effectType == EffectType.UTR_5_PRIME) //
				|| (effectType == EffectType.UTR_3_PRIME) //
				|| (effectType == EffectType.UTR_5_DELETED) //
				|| (effectType == EffectType.UTR_3_DELETED) //
		;
	}

	/**
	 * Create a list with only one element (this results)
	 * @return
	 */
	public List<ChangeEffectOld> newList() {
		List<ChangeEffectOld> list = new ArrayList<ChangeEffectOld>();
		list.add(clone());
		return list;
	}

	public void set(Marker marker, EffectType effectType, String message) {
		this.marker = marker;
		if (marker instanceof Exon) exon = (Exon) marker;
		this.effectType = effectType;
		this.message += message;
	}

	/**
	 * Set codon change. Calculate effect type based on codon changes (for SNPs ans MNPs)
	 * @param codonsOld
	 * @param codonsNew
	 * @param codonNum
	 * @param codonIndex
	 */
	public void setCodons(String codonsOld, String codonsNew, int codonNum, int codonIndex) {
		// Replace empty by "-"
		if (codonsOld.isEmpty()) codonsOld = "-";
		if (codonsNew.isEmpty()) codonsNew = "-";

		this.codonsOld = codonsOld;
		this.codonsNew = codonsNew;
		this.codonNum = codonNum;
		this.codonIndex = codonIndex;

		CodonTable codonTable = marker.codonTable();
		boolean indel = codonsOld.length() != codonsNew.length();

		// Calculate amino acids
		if (codonsOld.equals("-")) {
			aaOld = "-";
			indel = true;
		} else {
			aaOld = codonTable.aa(codonsOld);
			codonDegeneracy = codonTable.degenerate(codonsOld, codonIndex); // Calculate codon degeneracy
		}

		if (codonsNew.equals("-")) {
			aaNew = "-";
			indel = true;
		} else aaNew = codonTable.aa(codonsNew);

		if (!indel) {
			// SNM and MNP effects
			if (aaOld.equals(aaNew)) {
				// Same AA: Synonymous coding
				if ((codonNum == 0) && codonTable.isStartFirst(codonsOld)) {
					// It is in the first codon (which also is a start codon)
					if (codonTable.isStartFirst(codonsNew)) effectType = EffectType.SYNONYMOUS_START; // The new codon is also a start codon => SYNONYMOUS_START
					else effectType = EffectType.START_LOST; // The AA is the same, but the codon is not a start codon => start lost
				} else if (codonTable.isStop(codonsOld)) {
					// Stop codon
					if (codonTable.isStop(codonsNew)) effectType = EffectType.SYNONYMOUS_STOP; // New codon is also a stop => SYNONYMOUS_STOP
					else effectType = EffectType.STOP_LOST; // New codon is not a stop, the we've lost a stop
				} else {
					// All other cases are just SYNONYMOUS_CODING
					effectType = EffectType.SYNONYMOUS_CODING;
				}
			} else {
				// Different AA: Non-synonymous coding
				if ((codonNum == 0) && codonTable.isStartFirst(codonsOld)) {
					// It is in the first codon (which also is a start codon)
					if (codonTable.isStartFirst(codonsNew)) effectType = EffectType.NON_SYNONYMOUS_START; // Non-synonymous mutation on first codon => start lost
					else effectType = EffectType.START_LOST; // Non-synonymous mutation on first codon => start lost
				} else if (codonTable.isStop(codonsOld)) {
					// Stop codon
					if (codonTable.isStop(codonsNew)) effectType = EffectType.NON_SYNONYMOUS_STOP; // Notice: This should never happen! (for some reason I removed this comment at some point and that create some confusion): http://www.biostars.org/post/show/51352/in-snpeff-impact-what-is-difference-between-stop_gained-and-non-synonymous_stop/
					else effectType = EffectType.STOP_LOST;
				} else if (codonTable.isStop(codonsNew)) effectType = EffectType.STOP_GAINED;
				else {
					// All other cases are just NON_SYN
					effectType = EffectType.NON_SYNONYMOUS_CODING;
				}
			}
		} else {
			// Only change effect in some cases
			if ((codonNum == 0) && codonTable.isStartFirst(codonsOld) && !codonTable.isStartFirst(codonsNew)) effectType = EffectType.START_LOST;
			else if (codonTable.isStop(codonsOld) && !codonTable.isStop(codonsNew)) effectType = EffectType.STOP_LOST;
			else if (!codonTable.isStop(codonsOld) && codonTable.isStop(codonsNew)) effectType = EffectType.STOP_GAINED;
		}
	}

	/**
	 * Set values for codons around change.
	 * @param codonsLeft
	 * @param codonsRight
	 */
	public void setCodonsAround(String codonsLeft, String codonsRight) {
		codonsAroundOld = codonsLeft.toLowerCase() + codonsOld.toUpperCase() + codonsRight.toLowerCase();
		codonsAroundNew = codonsLeft.toLowerCase() + codonsNew.toUpperCase() + codonsRight.toLowerCase();

		// Amino acids surrounding the ones changed
		CodonTable codonTable = marker.codonTable();
		String aasLeft = codonTable.aa(codonsLeft);
		String aasRigt = codonTable.aa(codonsRight);
		aasAroundOld = aasLeft.toLowerCase() + aaOld.toUpperCase() + aasRigt.toLowerCase();
		aasAroundNew = aasLeft.toLowerCase() + aaNew.toUpperCase() + aasRigt.toLowerCase();
	}

	public void setExon(Exon exon) {
		this.exon = exon;
	}

	@Override
	public String toString() {
		// Get data to show
		String geneId = "", geneName = "", bioType = "", transcriptId = "", exonId = "", customId = "";
		int exonRank = -1;

		if (marker != null) {
			// Gene Id, name and biotype
			Gene gene = (Gene) marker.findParent(Gene.class);

			// CDS size info
			Transcript tr;
			if (exon != null) tr = (Transcript) exon.findParent(Transcript.class);
			else tr = (Transcript) marker.findParent(Transcript.class);

			if (gene != null) {
				geneId = gene.getId();
				geneName = gene.getGeneName();
				bioType = (tr != null ? tr.getBioType() : "");

				// Make one more effort to show whether gene is protein coding or not
				if ((bioType.isEmpty()) && gene.getGenome().hasCodingInfo()) bioType = (gene.isProteinCoding() ? "coding" : "non-coding");
			}

			// Update trId
			if (tr != null) transcriptId = tr.getId();

			// Exon rank information
			if (exon != null) {
				exonId = exon.getId();
				exonRank = exon.getRank();
			}

			// Regulation
			if (isRegulation()) {
				bioType = ((Regulation) marker).getCellType();
			}
		}

		// Add seqChage's ID
		if (!seqChange.getId().isEmpty()) customId += seqChange.getId();

		// Add custom markers
		if ((marker != null) && (marker instanceof Custom)) customId += (customId.isEmpty() ? "" : ";") + marker.getId();

		// CDS length
		int cdsSize = getCdsLength();

		String errWarn = error + (error.isEmpty() ? "" : "|") + warning;
		return errWarn //		
				+ "\t" + geneId //
				+ "\t" + geneName //
				+ "\t" + bioType //
				+ "\t" + transcriptId //
				+ "\t" + exonId //
				+ "\t" + (exonRank >= 0 ? exonRank : "") //
				+ "\t" + effect(false, false, false) //
				+ "\t" + ((aaOld.length() + aaNew.length()) > 0 ? aaOld + "/" + aaNew : "") //
				+ "\t" + ((codonsOld.length() + codonsNew.length()) > 0 ? codonsOld + "/" + codonsNew : "") //
				+ "\t" + (codonNum >= 0 ? (codonNum + 1) : "") //
				+ "\t" + (codonDegeneracy >= 0 ? codonDegeneracy + "" : "") //
				+ "\t" + (cdsSize >= 0 ? cdsSize : "") //
				+ "\t" + (codonsAroundOld.length() > 0 ? codonsAroundOld + " / " + codonsAroundNew : "") //
				+ "\t" + (aasAroundOld.length() > 0 ? aasAroundOld + " / " + aasAroundNew : "") //
				+ "\t" + customId //
		;
	}

	/**
	 * Get the simplest string describing the effect (this is mostly used for testcases)
	 * @param shortFormat
	 * @return
	 */
	public String toStringSimple(boolean shortFormat) {
		// Get data to show
		String transcriptId = "", exonId = "";

		if (marker != null) {
			Transcript tr = (Transcript) marker.findParent(Transcript.class);
			if (tr != null) transcriptId = tr.getId();

			if (exon != null) exonId = exon.getId();
		}

		String eff = effect(shortFormat, true, true);
		if (eff.length() > 0) return eff;
		if (exonId.length() > 0) return exonId;
		if (transcriptId.length() > 0) return transcriptId;
		return "NO EFFECT";
	}

}