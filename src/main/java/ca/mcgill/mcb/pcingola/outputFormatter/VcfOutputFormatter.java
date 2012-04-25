package ca.mcgill.mcb.pcingola.outputFormatter;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;

import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Regulation;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.FunctionalClass;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;

/**
 * Formats output as VCF
 * 
 * @author pcingola
 */
public class VcfOutputFormatter extends OutputFormatter {

	boolean needAddInfo = false;
	boolean needAddHeader = true;

	public VcfOutputFormatter() {
		super();
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
			vcfFile.addHeader(newHeaderLine);

		needAddHeader = false;
	}

	/**
	 * Add effects to INFO field
	 */
	protected void addInfo(VcfEntry vcfEntry) {
		// No effects to show?
		if (changeEffects.isEmpty()) return;

		//---
		// Calculate all effects and genes
		//---
		HashSet<String> effs = new HashSet<String>();
		for (ChangeEffect changeEffect : changeEffects) {
			// If it is not filtered out by changeEffectResutFilter  => Show it
			if ((changeEffectResutFilter == null) || (!changeEffectResutFilter.filter(changeEffect))) {
				StringBuilder effBuff = new StringBuilder();

				// Add effect
				effBuff.append(changeEffect.effect(true, false, false));
				effBuff.append("(");

				// Add effect impact
				effBuff.append(changeEffect.getEffectImpact());
				effBuff.append("|");

				// Add functional class
				FunctionalClass fc = changeEffect.getFunctionalClass();
				effBuff.append(fc == FunctionalClass.NONE ? "" : fc.toString()); // Show only if it is not empty
				effBuff.append("|");

				// Codon change
				effBuff.append(changeEffect.getCodonChange());
				effBuff.append("|");

				// Add amino acid change
				effBuff.append(changeEffect.getAaChangeHgvs());
				effBuff.append("|");

				// Add gene info
				Gene gene = changeEffect.getGene();
				if (gene != null) {
					// Protein coding gene?
					String coding = "";
					if (gene.getGenome().hasCodingInfo()) coding = (gene.isProteinCoding() ? "CODING" : "NON_CODING");

					effBuff.append(gene.getGeneName());
					effBuff.append("|");
					effBuff.append(gene.getBioType());
					effBuff.append("|");
					effBuff.append(coding);
					effBuff.append("|");
				} else if (changeEffect.isRegulation()) {
					Regulation reg = (Regulation) changeEffect.getMarker();
					effBuff.append("|" + reg.getCellType() + "||");
				} else effBuff.append("|||");

				// Add transcript info
				Transcript tr = changeEffect.getTranscript();
				if (tr != null) effBuff.append(tr.getId());
				effBuff.append("|");

				// Add exon info
				Exon ex = changeEffect.getExon();
				if (ex != null) effBuff.append(ex.getId());

				// Errors or warnings (this is the last thing in the list)
				if (!changeEffect.getWarning().isEmpty()) effBuff.append("|" + changeEffect.getWarning());
				if (!changeEffect.getError().isEmpty()) effBuff.append("|" + changeEffect.getError());

				effBuff.append(")");

				// Get effect
				effs.add(effBuff.toString());
			}
		}

		//---
		// Add effects (sorted)
		//---
		ArrayList<String> listEffs = new ArrayList<String>(effs);
		Collections.sort(listEffs);
		StringBuffer sbEffs = new StringBuffer();
		for (String eff : listEffs)
			sbEffs.append(eff + ",");

		if (sbEffs.length() > 0) sbEffs.deleteCharAt(sbEffs.length() - 1); // Remove last comma

		// Add 'EFF' info field
		vcfEntry.addInfo("EFF", sbEffs.toString());

		needAddInfo = false; // Don't add info twice
	}

	/**
	 * Finish up section
	 * @param marker
	 */
	@Override
	public String endSection(Marker marker) {
		// Ignore other markers (e.g. seqChanges)
		if (marker instanceof VcfEntry) return super.endSection(marker);
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
		newLines.add("##INFO=<ID=EFF,Number=.,Type=String,Description=\"Predicted effects for this variant.Format: 'Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_change | Gene_Name | Gene_BioType | Coding | Transcript | Exon [ | ERRORS | WARNINGS ] )' \">");
		return newLines;
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
		return vcfFile.getHeader();
	}
}
