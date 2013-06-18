package ca.mcgill.mcb.pcingola.snpEffect;

import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * A loss of function (LOF) entry in a VCF file 
 * 
 * @author pcingola
 */
public class LossOfFunctionEntry {

	String geneName, geneId;
	int transcriptInGene;
	double percentOfTranscriptsAffected;

	public LossOfFunctionEntry(Gene gene, double percentOfTranscriptsAffected) {
		geneName = gene.getGeneName();
		geneId = gene.getId();
		transcriptInGene = gene.numChilds();
		this.percentOfTranscriptsAffected = percentOfTranscriptsAffected;
	}

	public LossOfFunctionEntry(String vcfLofStr) {
		String fields[] = vcfLofStr.split("\\|");
		geneName = fields[0];
		geneId = fields[1];
		transcriptInGene = Gpr.parseIntSafe(fields[2]);
		percentOfTranscriptsAffected = Gpr.parseDoubleSafe(fields[3]);
	}

	public LossOfFunctionEntry(String geneName, String geneId, int transcriptInGene, double percentOfTranscriptsAffected) {
		this.geneName = geneName;
		this.geneId = geneId;
		this.transcriptInGene = transcriptInGene;
		this.percentOfTranscriptsAffected = percentOfTranscriptsAffected;
	}

	public String getGeneId() {
		return geneId;
	}

	public String getGeneName() {
		return geneName;
	}

	public double getPercentOfTranscriptsAffected() {
		return percentOfTranscriptsAffected;
	}

	public int getTranscriptInGene() {
		return transcriptInGene;
	}

	@Override
	public String toString() {
		return String.format("%s|%s|%d|%.2f", geneName, geneId, transcriptInGene, percentOfTranscriptsAffected);
	}
}
