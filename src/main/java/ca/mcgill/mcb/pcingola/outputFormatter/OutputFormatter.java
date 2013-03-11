package ca.mcgill.mcb.pcingola.outputFormatter;

import java.util.ArrayList;

import ca.mcgill.mcb.pcingola.filter.ChangeEffectFilter;
import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Genome;
import ca.mcgill.mcb.pcingola.interval.Intron;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;

/**
 * Formats output
 * How is this used:
 *    - newSection();   // Create a new 'section' on the output format (e.g. a new seqChange)
 *    - add();			// Add all changes related to this section (i.e. all changes related to this seqChange)
 *    - endSection();	// Output all changes related to this section (output header if needed), clean up list of changes
 *     
 * @author pcingola
 */
public abstract class OutputFormatter {

	boolean supressOutput; // Do not print anything
	boolean useHgvs; // Use HGVS notation
	boolean useGeneId; // Use Gene ID instead of gene name
	boolean useSequenceOntolgy; // Use Sequence Ontolgy terms
	boolean useOicr; // Use OICR tag
	boolean showHeader = true; // Show header information
	int sectionNum = 0;
	int outOffset = 1;
	String commandLineStr;
	String version;
	String chrStr;
	Marker section;
	ChangeEffectFilter changeEffectResutFilter = null; // Filter prediction results
	ArrayList<ChangeEffect> changeEffects;

	public static String idChain(Marker marker) {
		return idChain(marker, ";", true);
	}

	/**
	 * A list of all IDs and parent IDs until chromosome
	 * @param m
	 * @return
	 */
	public static String idChain(Marker marker, String separator, boolean useGeneId) {
		if (marker == null) return "";

		StringBuilder sb = new StringBuilder();

		for (Marker m = marker; (m != null) && !(m instanceof Chromosome) && !(m instanceof Genome); m = m.getParent()) {
			if (sb.length() > 0) sb.append(separator);
			switch (m.getType()) {
			case EXON:
				Transcript tr = (Transcript) m.getParent();
				sb.append("exon_" + ((Exon) m).getRank() + "_" + tr.numChilds());
				break;

			case INTRON:
				sb.append("intron_" + ((Intron) m).getRank());
				break;

			case GENE:
				Gene g = (Gene) m;
				sb.append(useGeneId ? m.getId() : g.getGeneName());
				sb.append(separator + g.getBioType());
				break;

			case TRANSCRIPT:
				sb.append(m.getId());
				sb.append(separator + ((Transcript) m).getBioType());
				break;

			case CHROMOSOME:
			case INTERGENIC:
				sb.append(m.getId());
				break;

			default:
				break;
			}
		}

		// Empty? Add ID
		if (sb.length() <= 0) sb.append(marker.getId());

		// Prepend type
		sb.insert(0, marker.getClass().getSimpleName() + separator);

		return sb.toString();
	}

	public OutputFormatter() {
		changeEffects = new ArrayList<ChangeEffect>();
	}

	/**
	 * Add effects to list
	 * @param changeEffects
	 */
	public void add(ChangeEffect changeEffect) {
		// Passes the filter? => Add
		if ((changeEffectResutFilter == null) || (!changeEffectResutFilter.filter(changeEffect))) changeEffects.add(changeEffect);
	}

	@Override
	public OutputFormatter clone() {
		OutputFormatter newOutputFormatter = null;
		try {
			// Create a new formatter. We cannot use the same output formatter for all workers
			newOutputFormatter = this.getClass().newInstance();
			newOutputFormatter.supressOutput = supressOutput;
			newOutputFormatter.showHeader = showHeader;
			newOutputFormatter.sectionNum = sectionNum;
			newOutputFormatter.outOffset = outOffset;
			newOutputFormatter.commandLineStr = commandLineStr;
			newOutputFormatter.version = version;
			newOutputFormatter.chrStr = chrStr;
			newOutputFormatter.section = section;
			newOutputFormatter.changeEffectResutFilter = changeEffectResutFilter;
		} catch (Exception e) {
			throw new RuntimeException(e);
		}

		return newOutputFormatter;
	}

	/**
	 * Finish up section
	 * @param marker
	 */
	public String endSection(Marker marker) {
		StringBuilder sb = null;

		if (!supressOutput) {
			sb = new StringBuilder();
			// Add header?
			if (showHeader && (sectionNum == 0)) {
				String header = toStringHeader();
				if (!header.isEmpty()) {
					sb.append(header);
					sb.append("\n");
				}
			}

			// Add current line
			sb.append(toString());
		}

		sectionNum++;
		changeEffects.clear();

		if (supressOutput) return null;
		return sb.toString();
	}

	/**
	 * End this section and print results
	 * @param marker
	 */
	public void printSection(Marker marker) {
		String outStr = endSection(marker);
		if ((outStr != null) && (!outStr.isEmpty())) System.out.println(outStr);
	}

	public void setChangeEffectResutFilter(ChangeEffectFilter changeEffectResutFilter) {
		this.changeEffectResutFilter = changeEffectResutFilter;
	}

	public void setChrStr(String chrStr) {
		this.chrStr = chrStr;
	}

	public void setCommandLineStr(String commandLineStr) {
		this.commandLineStr = commandLineStr;
	}

	public void setOutOffset(int outOffset) {
		this.outOffset = outOffset;
	}

	public void setShowHeader(boolean showHeader) {
		this.showHeader = showHeader;
	}

	public void setSupressOutput(boolean supressOutput) {
		this.supressOutput = supressOutput;
	}

	public void setUseGeneId(boolean useGeneId) {
		this.useGeneId = useGeneId;
	}

	public void setUseHgvs(boolean useHgvs) {
		this.useHgvs = useHgvs;
	}

	public void setUseOicr(boolean useOicr) {
		this.useOicr = useOicr;
	}

	public void setUseSequenceOntolgy(boolean useSequenceOntolgy) {
		this.useSequenceOntolgy = useSequenceOntolgy;
	}

	public void setVersion(String version) {
		this.version = version;
	}

	/**
	 * Starts a new section
	 * @param section
	 */
	public void startSection(Marker marker) {
		section = marker;
	}

	@Override
	public String toString() {
		throw new RuntimeException("Method toString() must be overridden!");
	}

	/**
	 * Show header
	 */
	protected abstract String toStringHeader();
}
