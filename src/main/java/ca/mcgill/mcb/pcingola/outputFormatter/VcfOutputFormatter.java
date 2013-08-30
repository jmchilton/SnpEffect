package ca.mcgill.mcb.pcingola.outputFormatter;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;

import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Genome;
import ca.mcgill.mcb.pcingola.interval.Intron;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Regulation;
import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.FunctionalClass;
import ca.mcgill.mcb.pcingola.snpEffect.LossOfFunction;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.vcf.VcfEffect;
import ca.mcgill.mcb.pcingola.vcf.VcfEffect.FormatVersion;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;

/**
 * Formats output as VCF
 * 
 * @author pcingola
 */
public class VcfOutputFormatter extends OutputFormatter {

	public static final boolean debug = false;

	//	public static final String GATK_ACCEPTED_VERSION = "2.0.5"; // GATK refuses to run if we report another version, so we have to lie...
	public static final String VCF_INFO_OICR_NAME = "OICR";

	boolean needAddInfo = false;
	boolean needAddHeader = true;
	FormatVersion formatVersion = VcfEffect.FormatVersion.FORMAT_SNPEFF_4;
	List<VcfEntry> vcfEntries;
	boolean lossOfFunction;
	boolean gatk;
	Genome genome;

	/**
	 * This constructor is used mostly by 'clone()' method
	 */
	protected VcfOutputFormatter() {
		super();
	}

	public VcfOutputFormatter(Genome genome) {
		super();
		this.genome = genome;
	}

	/**
	 * Add all vcf entries to a list (used only for debugging and test-cases)
	 * @param vcfEntries
	 */
	public VcfOutputFormatter(List<VcfEntry> vcfEntries) {
		super();
		this.vcfEntries = vcfEntries;
	}

	/**
	 * Add header
	 * @param vcfEntry
	 */
	protected void addHeader() {
		VcfEntry vcfEntry = (VcfEntry) section;

		// Get header
		VcfFileIterator vcfFile = vcfEntry.getVcfFileIterator();

		// Add new lines
		for (String newHeaderLine : getNewHeaderLines())
			vcfFile.getVcfHeader().addLine(newHeaderLine);

		needAddHeader = false;
	}

	/**
	 * Add effects to INFO field
	 */
	protected void addInfo(VcfEntry vcfEntry) {
		// No effects to show?
		if (changeEffects.isEmpty()) return;

		// Do all effects have warnings or errors (this could produce an empty 'EFF' field in GATK mode?
		boolean allWarnings = false;
		if (gatk) {
			allWarnings = changeEffects.size() > 0;
			for (ChangeEffect changeEffect : changeEffects)
				allWarnings &= (changeEffect.hasError() || changeEffect.hasWarning());
		}

		//---
		// Calculate all effects and genes
		//---
		HashSet<String> effs = new HashSet<String>();
		HashSet<String> oicr = (useOicr ? new HashSet<String>() : null);
		for (ChangeEffect changeEffect : changeEffects) {
			// In GATK mode, skip changeEffects having errors or warnings (unless ALL effects have warnings)
			if (gatk && !allWarnings && (changeEffect.hasError() || changeEffect.hasWarning())) continue;

			// If it is not filtered out by changeEffectResutFilter  => Show it
			if ((changeEffectResutFilter == null) || (!changeEffectResutFilter.filter(changeEffect))) {
				StringBuilder effBuff = new StringBuilder();

				// Add effect
				effBuff.append(changeEffect.effect(true, false, false, useSequenceOntolgy));
				effBuff.append("(");

				// Add effect impact
				effBuff.append(changeEffect.getEffectImpact());
				effBuff.append("|");

				// Add functional class
				FunctionalClass fc = changeEffect.getFunctionalClass();
				effBuff.append(fc == FunctionalClass.NONE ? "" : fc.toString()); // Show only if it is not empty
				effBuff.append("|");

				// Codon change
				String codonChange = changeEffect.getCodonChange();
				if (!codonChange.isEmpty()) effBuff.append(codonChange);
				else if (changeEffect.getDistance() >= 0) effBuff.append(changeEffect.getDistance());
				effBuff.append("|");

				// Add HGVS (amino acid change) 
				if (useHgvs) effBuff.append(changeEffect.getHgvs());
				else effBuff.append(changeEffect.getAaChangeHgsv());
				effBuff.append("|");

				// Add amino acid length
				if (formatVersion != FormatVersion.FORMAT_SNPEFF_2) { // This field is not in format version 2
					int aalen = changeEffect.getAaLength();
					effBuff.append(aalen >= 0 ? aalen : "");
					effBuff.append("|");
				}

				// Add gene info
				Gene gene = changeEffect.getGene();
				Transcript tr = changeEffect.getTranscript();
				if (gene != null) {
					// Gene name
					effBuff.append(vcfInfoSafeString(useGeneId ? gene.getId() : gene.getGeneName()));
					effBuff.append("|");

					// Transcript ID
					effBuff.append(tr != null ? tr.getBioType() : "");
					effBuff.append("|");

					// Protein coding gene?
					String coding = "";
					if (gene.getGenome().hasCodingInfo()) coding = (gene.isProteinCoding() ? ChangeEffect.Coding.CODING.toString() : ChangeEffect.Coding.NON_CODING.toString());
					effBuff.append(coding);
					effBuff.append("|");
				} else if (changeEffect.isRegulation()) {
					Regulation reg = (Regulation) changeEffect.getMarker();
					effBuff.append("|" + reg.getCellType() + "||");
				} else if (changeEffect.isCustom()) {
					Marker m = changeEffect.getMarker();
					if (m != null) effBuff.append("|" + m.getId() + "||");
					else effBuff.append("|||");
				} else effBuff.append("|||");

				// Add transcript info
				if (tr != null) effBuff.append(vcfInfoSafeString(tr.getId()));
				effBuff.append("|");

				// Add exon (or intron) rank info
				Exon ex = changeEffect.getExon();
				int rank = -1;
				if (ex != null) rank = ex.getRank();
				else {
					// Do we have an intron?
					Intron intron = changeEffect.getIntron();
					if (intron != null) rank = intron.getRank();
				}
				effBuff.append(rank >= 0 ? rank : "");

				// Add genotype (or genotype difference) for this effect
				if (formatVersion == FormatVersion.FORMAT_SNPEFF_4) {
					effBuff.append("|");
					effBuff.append(changeEffect.getGenotype());
				}

				//---
				// Errors or warnings (this is the last thing in the list)
				//---
				if (changeEffect.hasError() || changeEffect.hasWarning()) {
					StringBuilder err = new StringBuilder();

					// Add warnings
					if (!changeEffect.getWarning().isEmpty()) err.append(changeEffect.getWarning());

					// Add errors
					if (!changeEffect.getError().isEmpty()) {
						if (err.length() > 0) err.append("+");
						err.append(changeEffect.getError());
					}

					effBuff.append("|");
					effBuff.append(err);
				}
				effBuff.append(")");

				//---
				// Add effect
				//---
				if (!effs.add(effBuff.toString())) {
					if (debug) {
						// Effect has already been added? Something is wrong, the information should be unique for each effect
						StringBuilder sb = new StringBuilder();
						sb.append("--------------------------------------------------------------------------------\n");
						sb.append("VCF Entry   :\t" + vcfEntry + "\n");
						sb.append("REPEAT (VCF):\t" + effBuff + "\n");
						sb.append("REPEAT (TXT):\t" + changeEffect + "\n");
						sb.append("All    (VCF):\n");
						for (String ce : effs)
							sb.append("\t" + ce + "\n");
						sb.append("All    (TXT):\n");
						for (ChangeEffect ce : changeEffects)
							sb.append("\t" + ce + "\n");
						sb.append("--------------------------------------------------------------------------------\n");
						Gpr.debug("WARNING: Repeated effect!\n" + sb);
					}
				}

				//---
				// Add OICR data
				// 
				//---
				if (useOicr && (tr != null)) {
					StringBuilder sb = new StringBuilder();
					SeqChange seqChange = changeEffect.getSeqChange();

					// Get cDNA position
					int pos = tr.isStrandMinus() ? seqChange.getStart() : seqChange.getEnd(); // First base in cDNA
					int cdnaIdx = tr.cDnaBaseNumber(pos); // Which cDNA base number?
					if (cdnaIdx > 0) sb.append("(" + tr.getId() + "|" + cdnaIdx + ")");

					oicr.add(sb.toString());
				}
			}
		}

		//---
		// Add data to INFO fields
		//---

		// Add 'EFF' info field
		String effStr = toStringVcfInfo(effs);
		if (!effStr.isEmpty()) vcfEntry.addInfo(VcfEffect.VCF_INFO_EFF_NAME, effStr);

		// Add 'OICR' info field
		if (useOicr && (oicr.size() > 0)) {
			String oicrInfo = toStringVcfInfo(oicr);
			if (!oicrInfo.isEmpty()) vcfEntry.addInfo(VCF_INFO_OICR_NAME, oicrInfo);
		}

		// Add LOF info?
		if (lossOfFunction) {
			// Perform LOF analysis and add annotations
			LossOfFunction lof = new LossOfFunction(config, changeEffects);
			if (lof.isLof()) vcfEntry.addInfo(LossOfFunction.VCF_INFO_LOF_NAME, lof.toStringVcfLof());
			if (lof.isNmd()) vcfEntry.addInfo(LossOfFunction.VCF_INFO_NMD_NAME, lof.toStringVcfNmd());
		}

		needAddInfo = false; // Don't add info twice
	}

	@Override
	public OutputFormatter clone() {
		try {
			VcfOutputFormatter newOutputFormatter = (VcfOutputFormatter) super.clone();
			newOutputFormatter.formatVersion = formatVersion;
			newOutputFormatter.needAddInfo = needAddInfo;
			newOutputFormatter.needAddHeader = needAddHeader;
			newOutputFormatter.lossOfFunction = lossOfFunction;
			newOutputFormatter.genome = genome;
			newOutputFormatter.config = config;
			return newOutputFormatter;
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * Finish up section
	 * @param marker
	 */
	@Override
	public String endSection(Marker marker) {
		// Ignore other markers (e.g. seqChanges)
		if (marker instanceof VcfEntry) {
			if (vcfEntries != null) vcfEntries.add((VcfEntry) marker);
			return super.endSection(marker);
		}
		return null;
	}

	/**
	 * New lines to be added to header
	 * @return
	 */
	public List<String> getNewHeaderLines() {
		ArrayList<String> newLines = new ArrayList<String>();

		newLines.add("##SnpEffVersion=\"" + version + "\"");
		newLines.add("##SnpEffCmd=\"" + commandLineStr + "\"");

		// Fields changed in different format versions
		if (formatVersion == FormatVersion.FORMAT_SNPEFF_2) newLines.add("##INFO=<ID=EFF,Number=.,Type=String,Description=\"Predicted effects for this variant.Format: 'Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_change| Gene_Name | Transcript_BioType | Gene_Coding | Transcript_ID | Exon [ | ERRORS | WARNINGS ] )' \">");
		else if (formatVersion == FormatVersion.FORMAT_SNPEFF_3) newLines.add("##INFO=<ID=EFF,Number=.,Type=String,Description=\"Predicted effects for this variant.Format: 'Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_change| Amino_Acid_length | Gene_Name | Transcript_BioType | Gene_Coding | Transcript_ID | Exon [ | ERRORS | WARNINGS ] )' \">");
		else newLines.add("##INFO=<ID=EFF,Number=.,Type=String,Description=\"Predicted effects for this variant.Format: 'Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_change| Amino_Acid_length | Gene_Name | Transcript_BioType | Gene_Coding | Transcript_ID | Exon  | GenotypeNum [ | ERRORS | WARNINGS ] )' \">");

		if (lossOfFunction) {
			newLines.add("##INFO=<ID=LOF,Number=.,Type=String,Description=\"Predicted loss of function effects for this variant. Format: 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected' \">");
			newLines.add("##INFO=<ID=NMD,Number=.,Type=String,Description=\"Predicted nonsense mediated decay effects for this variant. Format: 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected' \">");
		}

		if (useOicr) newLines.add("##INFO=<ID=OICR,Number=.,Type=String,Description=\"Format: ( Transcript | Distance from begining cDNA )\">");

		return newLines;
	}

	public void setGatk(boolean gatk) {
		this.gatk = gatk;
		if (gatk) formatVersion = VcfEffect.FormatVersion.FORMAT_SNPEFF_2;
	}

	public void setLossOfFunction(boolean lossOfFunction) {
		this.lossOfFunction = lossOfFunction;
	}

	@Override
	public void startSection(Marker marker) {
		// Ignore other markers (e.g. seqChanges)
		if (marker instanceof VcfEntry) super.startSection(marker);
		needAddInfo = true;
	}

	@Override
	public String toString() {
		VcfEntry vcfEntry = (VcfEntry) section;
		if (needAddInfo) addInfo(vcfEntry);
		return vcfEntry.toString();
	}

	/**
	 * Show header
	 */
	@Override
	protected String toStringHeader() {
		if (needAddHeader) addHeader(); // Add header lines

		VcfEntry vcfEntry = (VcfEntry) section;
		VcfFileIterator vcfFile = vcfEntry.getVcfFileIterator();
		return vcfFile.getVcfHeader().toString();
	}

	/**
	 * Convert a collection to a string usable in a VCF INFO field
	 * @param strs
	 * @return
	 */
	String toStringVcfInfo(Collection<String> strs) {
		// Sort strings
		ArrayList<String> list = new ArrayList<String>(strs);
		Collections.sort(list);

		// Add the all
		StringBuffer sb = new StringBuffer();
		for (String str : list)
			if (!str.isEmpty()) sb.append(str + ",");

		if (sb.length() > 0) sb.deleteCharAt(sb.length() - 1); // Remove last comma
		return sb.toString();
	}

	/**
	 * Create a string that is safe (i.e. valid) to add in an INFO field
	 */
	public String vcfInfoSafeString(String value) {
		if (value == null) return value;
		value = value.replaceAll("[ ,;|=()]", "_");
		return value;
	}
}
