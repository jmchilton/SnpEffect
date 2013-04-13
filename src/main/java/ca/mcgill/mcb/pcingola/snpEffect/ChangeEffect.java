package ca.mcgill.mcb.pcingola.snpEffect;

import java.util.ArrayList;
import java.util.List;

import ca.mcgill.mcb.pcingola.codons.CodonTable;
import ca.mcgill.mcb.pcingola.interval.Custom;
import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Intron;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.NextProt;
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
public class ChangeEffect implements Cloneable, Comparable<ChangeEffect> {

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
		, CUSTOM //
		, NEXT_PROT //
		;

		public String toSequenceOntology() {
			switch (this) {
			case CHROMOSOME:
				return "chromosome";
			case INTERGENIC:
				return "intergenic_region";
			case UPSTREAM:
				return "upstream_gene_variant";
			case UTR_5_PRIME:
				return "5_prime_UTR_variant";
			case UTR_5_DELETED:
				return "five_prime_UTR";
			case SPLICE_SITE_ACCEPTOR:
				return "splice_region_variant";
			case SPLICE_SITE_BRANCH:
				return "splice_region_variant";
			case SPLICE_SITE_BRANCH_U12:
				return "splice_region_variant";
			case SPLICE_SITE_DONOR:
				return "splice_region_variant";
			case START_LOST:
				return "initiator_codon_variant";
			case SYNONYMOUS_START:
				return "initiator_codon_variant";
			case NON_SYNONYMOUS_START:
				return "initiator_codon_variant";
			case TRANSCRIPT:
				return "nc_transcript_variant";
			case EXON:
				return "non_coding_exon_variant";
			case EXON_DELETED:
				return "exon_lost";
			case NON_SYNONYMOUS_CODING:
				return "missense";
			case SYNONYMOUS_CODING:
				return "synonymous_variant";
			case FRAME_SHIFT:
				return "frameshift_variant";
			case CODON_CHANGE:
				return "coding_sequence_variant";
			case CODON_INSERTION:
				return "inframe_insertion";
			case CODON_CHANGE_PLUS_CODON_INSERTION:
				return "inframe_insertion";
			case CODON_DELETION:
				return "inframe_deletion";
			case CODON_CHANGE_PLUS_CODON_DELETION:
				return "inframe_deletion";
			case STOP_GAINED:
				return "stop_gained";
			case SYNONYMOUS_STOP:
				return "stop_retained_variant";
			case NON_SYNONYMOUS_STOP:
				return "stop_retained_variant";
			case STOP_LOST:
				return "stop_lost";
			case INTRON:
				return "intron_variant";
			case UTR_3_PRIME:
				return "3_prime_UTR_variant";
			case UTR_3_DELETED:
				return "";
			case DOWNSTREAM:
				return "downstream_gene_variant";
			case INTRON_CONSERVED:
				return "intron_variant";
			case INTERGENIC_CONSERVED:
				return "intergenic_variant";
			case INTRAGENIC:
				return "intergenic_variant";
			case REGULATION:
				return "regulatory_region_variant";
			case RARE_AMINO_ACID:
				return "non_conservative_missense_variant";

			case START_GAINED:
			case MICRO_RNA:
			case NONE:
			case GENE:
			case CDS:
			case GENOME:
			case CUSTOM:
				return this.toString().toLowerCase(); // Just a wild guess ... this should probably throw an Exception

			default:
				throw new RuntimeException("Sequence Ontology term not found for EffectType '" + this + "'");
			}
		}
	};

	/**
	 * Errors for change effect
	 * @author pcingola
	 *
	 */
	public enum ErrorType {
		ERROR_CHROMOSOME_NOT_FOUND //
		, ERROR_OUT_OF_CHROMOSOME_RANGE //
		, ERROR_OUT_OF_EXON //
		, ERROR_MISSING_CDS_SEQUENCE //
	}

	/**
	 * This class is only getFused for SNPs
	 */
	public enum FunctionalClass {
		NONE, SILENT, MISSENSE, NONSENSE
	}

	public enum WarningType {
		WARNING_SEQUENCE_NOT_AVAILABLE //
		, WARNING_REF_DOES_NOT_MATCH_GENOME //
	};

	static final boolean COMPATIBLE_v1_8 = true; // Activate this in order to get the same out as version 1.8. This is only for testing & debugging 

	SeqChange seqChange = null;
	SeqChange seqChangeRef = null;
	EffectType effectType = EffectType.NONE;
	EffectImpact effectImpact = null;
	Marker marker = null;
	String error = "", warning = "", message = ""; // Any message, warning or error?
	String codonsOld = "", codonsNew = ""; // Codon change information
	String codonsAroundOld = "", codonsAroundNew = ""; // Codons around
	int codonNum = -1; // Codon number (negative number mens 'information not available')
	int codonIndex = -1; // Index within a codon (negative number mens 'information not available')
	int codonDegeneracy = -1; // Codon degeneracy (negative number mens 'information not available')
	String aaOld = "", aaNew = ""; // Amino acid changes
	String aasAroundOld = "", aasAroundNew = ""; // Amino acids around

	/**
	 *  An empty list of results;
	 * @return
	 */
	public static List<ChangeEffect> emptyResults() {
		return new ArrayList<ChangeEffect>();
	}

	public ChangeEffect(SeqChange seqChange) {
		this.seqChange = seqChange;
	}

	public ChangeEffect(SeqChange seqChange, SeqChange seqChangeRef) {
		this.seqChange = seqChange;
		this.seqChangeRef = seqChangeRef;
	}

	public void addError(ErrorType err) {
		error += (error.isEmpty() ? "" : "+") + err;
	}

	public void addWarning(WarningType warn) {
		warning += (warning.isEmpty() ? "" : "+") + warn;
	}

	@Override
	public ChangeEffect clone() {
		try {
			return (ChangeEffect) super.clone();
		} catch (CloneNotSupportedException e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * Create a string for codon effect
	 * @param showAaChange : If true, include codon change, biotype, etc.
	 * @return
	 */
	String codonEffect(boolean showAaChange, boolean showBioType, boolean useSeqOntology) {
		String codonEffect = "";
		if ((marker == null) || (codonNum < 0)) return codonEffect;

		// Add codon effect
		codonEffect += getEffectTypeString(useSeqOntology);

		// Append codon change
		if (showAaChange) codonEffect += "(" + getAaChange() + ")";

		return codonEffect;
	}

	@Override
	public int compareTo(ChangeEffect changeEffect) {
		int comp = getEffectImpact().compareTo(changeEffect.getEffectImpact());
		if (comp != 0) return comp;

		comp = getEffectType().compareTo(changeEffect.getEffectType());
		if (comp != 0) return comp;

		if ((getMarker() != null) && (changeEffect.getMarker() != null)) return getMarker().compareTo(changeEffect.getMarker());

		return seqChange.compareTo(changeEffect.getSeqChange());
	}

	/**
	 * Show a string with overall effect
	 * @param shortFormat
	 * @param showAaChange
	 * @return
	 */
	public String effect(boolean shortFormat, boolean showAaChange, boolean showBioType, boolean useSeqOntology) {
		String e = "";
		String codonEffect = codonEffect(showAaChange, showBioType, useSeqOntology); // Codon effect

		// Create effect string
		if (!codonEffect.isEmpty()) e = codonEffect;
		else if (isRegulation()) e = getEffectTypeString(useSeqOntology) + "[" + ((Regulation) marker).getName() + "]";
		else if (isIntergenic() || isIntron() || isSpliceSite()) e = getEffectTypeString(useSeqOntology);
		else if (!message.isEmpty()) e = getEffectTypeString(useSeqOntology) + ": " + message;
		else if (marker == null) e = getEffectTypeString(useSeqOntology); // There are cases when no marker is associated (e.g. "Out of chromosome", "No such chromosome", etc.)
		else e = getEffectTypeString(useSeqOntology) + ": " + marker.getId();

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
	 * @return
	 */
	public String getAaChangeHgsv() {
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
	 * Get biotype
	 * @return
	 */
	public String getBiotype() {
		Gene gene = getGene();
		if (gene == null) return "";

		Transcript tr = getTranscript();
		if (tr != null) return tr.getBioType();
		else if (gene.getGenome().hasCodingInfo()) return (gene.isProteinCoding() ? "coding" : "non-coding");

		return "";
	}

	/**
	 * CDS length (negative if there is none)
	 * @return
	 */
	public int getCdsLength() {
		// CDS size info
		Transcript tr = getTranscript();
		if ((tr != null) && tr.isProteinCoding()) return tr.cds().length();
		return -1;
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

	public int getCodonNum() {
		return codonNum;
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

			case NEXT_PROT:
				if (marker == null) effectImpact = EffectImpact.MODIFIER;
				else if (((NextProt) marker).isHighlyConservedAaSequence()) effectImpact = EffectImpact.HIGH;
				else effectImpact = EffectImpact.LOW;
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

	/**
	 * Get Effect Type as a string
	 * @param useSeqOntology
	 * @return
	 */
	public String getEffectTypeString(boolean useSeqOntology) {
		if (effectType == null) return "";
		if (useSeqOntology) return effectType.toSequenceOntology();
		return effectType.toString();
	}

	public String getError() {
		return error;
	}

	/**
	 * Get exon (if any)
	 * @return
	 */
	public Exon getExon() {
		if (marker != null) {
			if (marker instanceof Exon) return (Exon) marker;
			return (Exon) marker.findParent(Exon.class);
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
		if (marker != null) {
			if (marker instanceof Gene) return (Gene) marker;
			return (Gene) marker.findParent(Gene.class);
		}
		return null;
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
		case NEXT_PROT:
		case TRANSCRIPT:
			if (isExon()) return EffectType.EXON.toString();
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

	/**
	 * Get genotype string
	 * @return
	 */
	public String getGenotype() {
		if (seqChange == null) return "";
		if (seqChangeRef != null) return seqChange.getGenotype() + "-" + seqChangeRef.getGenotype();
		return seqChange.getGenotype();
	}

	/**
	 * Change in HGVS notation
	 * References: http://www.hgvs.org/mutnomen/recs.html
	 * 
	 * @return
	 */
	public String getHgvs() {
		if (aaOld.isEmpty() && aaNew.isEmpty()) {
			if (codonNum >= 0) return "" + (codonNum + 1);
			return getHgvsNonCoding();
		}

		return getHgvsCoding();
	}

	/**
	 * Coding change in HGVS notation (amino acid changes)
	 * References: http://www.hgvs.org/mutnomen/recs.html
	 * 
	 * @return
	 */
	protected String getHgvsCoding() {
		// Codon numbering
		// HGVS: the translation initiator Methionine is numbered as +1
		int pos = codonNum + 1;

		// Synonymous changes
		if (aaOld.equals(aaNew)) {
			// HGVS: Description of so called "silent" changes in the format p.Leu54Leu (or p.L54L) is not allowed; descriptions 
			// 		 should be given at DNA level, it is non-informative and not unequivocal (there are five possibilities 
			// 		 at DNA level which may underlie p.Leu54Leu);  correct description has the format c.162C>G.
			pos = codonNum * 3 + codonIndex;
			return "c." + pos + seqChange.getReference() + ">" + seqChange.getChange();
		}

		// Convert to 3 letter code
		// HGVS: the three-letter amino acid code is prefered (see Discussion), with "*" designating a translation 
		// 		 termination codon; for clarity we this page describes changes using the three-letter amino acid
		String aaNew3 = aaNew;
		String aaOld3 = aaOld;
		if (marker != null) {
			CodonTable codonTable = marker.codonTable();
			aaNew3 = codonTable.aaThreeLetterCode(aaNew);
			aaOld3 = codonTable.aaThreeLetterCode(aaOld);

			// Start codon?
			if ((pos == 1) && codonTable.isStartFirst(aaOld)) {
				// HGVS: Currently, variants in the translation initiating Methionine (M1) are usually described as a substitution, 
				// 		 e.g. p.Met1Val. 
				//		 This is not correct. Either no protein is produced (p.0) or a new translation initiation site up- or downstream 
				//		 is used (e.g. p.Met1ValextMet-12 or p.Met1_Lys45del resp.). Unless experimental proof is available, it is probably 
				//		 best to report the effect on protein level as "p.Met1?" (unknown).
				return "p.Met1?";
			}
		}

		return "p." + aaOld3 + pos + aaNew3;
	}

	/**
	 * Coding change in HGVS notation (DNA changes)
	 * References: http://www.hgvs.org/mutnomen/recs.html
	 * 
	 * @return
	 */
	protected String getHgvsNonCoding() {
		return "";
	}

	/**
	 * Get intron (if any)
	 * @return
	 */
	public Intron getIntron() {
		if (marker != null) {
			if (marker instanceof Intron) return (Intron) marker;
			return (Intron) marker.findParent(Intron.class);
		}
		return null;
	}

	public Marker getMarker() {
		return marker;
	}

	public SeqChange getSeqChange() {
		return seqChange;
	}

	public Transcript getTranscript() {
		if (marker != null) {
			if (marker instanceof Transcript) return (Transcript) marker;
			return (Transcript) marker.findParent(Transcript.class);
		}
		return null;
	}

	public String getWarning() {
		return warning;
	}

	public boolean hasError() {
		return (error != null) && (!error.isEmpty());
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
		return (marker instanceof Exon) || (effectType == EffectType.EXON_DELETED);
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

	public boolean isNextProt() {
		return (effectType == EffectType.NEXT_PROT);
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
	public List<ChangeEffect> newList() {
		List<ChangeEffect> list = new ArrayList<ChangeEffect>();
		list.add(clone());
		return list;
	}

	public void set(Marker marker, EffectType effectType, String message) {
		this.marker = marker;
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

	public void setMarker(Marker marker) {
		this.marker = marker;
	}

	@Override
	public String toString() {
		return toString(false, false);
	}

	public String toString(boolean useSeqOntology, boolean useHgvs) {
		// Get data to show
		String geneId = "", geneName = "", bioType = "", transcriptId = "", exonId = "", customId = "";
		int exonRank = -1;

		if (marker != null) {
			// Gene Id, name and biotype
			Gene gene = getGene();
			Transcript tr = getTranscript();

			// CDS size info
			if (gene != null) {
				geneId = gene.getId();
				geneName = gene.getGeneName();
				bioType = getBiotype();
			}

			// Update trId
			if (tr != null) transcriptId = tr.getId();

			// Exon rank information
			Exon exon = getExon();
			if (exon != null) {
				exonId = exon.getId();
				exonRank = exon.getRank();
			}

			// Regulation
			if (isRegulation()) bioType = ((Regulation) marker).getCellType();
		}

		// Add seqChage's ID
		if (!seqChange.getId().isEmpty()) customId += seqChange.getId();

		// Add custom markers
		if ((marker != null) && (marker instanceof Custom)) customId += (customId.isEmpty() ? "" : ";") + marker.getId();

		// CDS length
		int cdsSize = getCdsLength();

		String errWarn = error + (error.isEmpty() ? "" : "|") + warning;

		String aaChange = "";
		if (useHgvs) aaChange = getHgvs();
		else aaChange = ((codonsOld.length() + codonsNew.length()) > 0 ? codonsOld + "/" + codonsNew : "");

		return errWarn //		
				+ "\t" + geneId //
				+ "\t" + geneName //
				+ "\t" + bioType //
				+ "\t" + transcriptId //
				+ "\t" + exonId //
				+ "\t" + (exonRank >= 0 ? exonRank : "") //
				+ "\t" + effect(false, false, false, useSeqOntology) //
				+ "\t" + aaChange //
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
		String transcriptId = "";
		Transcript tr = getTranscript();
		if (tr != null) transcriptId = tr.getId();

		String exonId = "";
		Exon exon = getExon();
		if (exon != null) exonId = exon.getId();

		String eff = effect(shortFormat, true, true, false);
		if (!eff.isEmpty()) return eff;
		if (!exonId.isEmpty()) return exonId;
		if (!transcriptId.isEmpty()) return transcriptId;

		return "NO EFFECT";
	}
}
