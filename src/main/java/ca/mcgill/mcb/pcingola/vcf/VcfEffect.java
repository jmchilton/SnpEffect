package ca.mcgill.mcb.pcingola.vcf;

import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * An 'EFF' entry in a vcf line
 * 
 * @author pablocingolani
 */
public class VcfEffect {

	ChangeEffect.EffectType effect;
	ChangeEffect.EffectImpact impact;
	ChangeEffect.FunctionalClass funClass;
	String codon;
	String aa;
	int aaLen;
	String gene;
	String bioType;
	ChangeEffect.Coding coding;
	String transcriptId;
	String exonId;

	public VcfEffect(String effStr) {
		parse(effStr);
	}

	public String getAa() {
		return aa;
	}

	public int getAaLen() {
		return aaLen;
	}

	public String getBioType() {
		return bioType;
	}

	public ChangeEffect.Coding getCoding() {
		return coding;
	}

	public String getCodon() {
		return codon;
	}

	public ChangeEffect.EffectType getEffect() {
		return effect;
	}

	public String getExonId() {
		return exonId;
	}

	public ChangeEffect.FunctionalClass getFunClass() {
		return funClass;
	}

	public String getGene() {
		return gene;
	}

	public ChangeEffect.EffectImpact getImpact() {
		return impact;
	}

	public String getTranscriptId() {
		return transcriptId;
	}

	void parse(String eff) {
		eff = eff.replace('(', ' '); // Replace all chars by spaces
		eff = eff.replace('|', ' ');
		eff = eff.replace(')', ' ');
		String effs[] = eff.split("\\s");

		// Parse each sub field
		effect = ChangeEffect.EffectType.valueOf(effs[0]);
		if (!effs[1].isEmpty()) impact = ChangeEffect.EffectImpact.valueOf(effs[1]);
		if (!effs[2].isEmpty()) funClass = ChangeEffect.FunctionalClass.valueOf(effs[2]);
		if (!effs[3].isEmpty()) codon = effs[3];
		if (!effs[4].isEmpty()) aa = effs[4];

		if (!effs[5].isEmpty()) aaLen = Gpr.parseIntSafe(effs[5]);
		else aaLen = 0;

		if ((effs.length >= 7) && !effs[6].isEmpty()) gene = effs[6];
		if ((effs.length >= 8) && !effs[7].isEmpty()) bioType = effs[7];
		if ((effs.length >= 9) && !effs[8].isEmpty()) coding = ChangeEffect.Coding.valueOf(effs[8]);
		if ((effs.length >= 10) && !effs[9].isEmpty()) transcriptId = effs[9];
		if ((effs.length >= 11) && !effs[10].isEmpty()) exonId = effs[10];
	}

	public void setAa(String aa) {
		this.aa = aa;
	}

	public void setAaLen(int aaLen) {
		this.aaLen = aaLen;
	}

	public void setBioType(String bioType) {
		this.bioType = bioType;
	}

	public void setCoding(ChangeEffect.Coding coding) {
		this.coding = coding;
	}

	public void setCodon(String codon) {
		this.codon = codon;
	}

	public void setEffect(ChangeEffect.EffectType effect) {
		this.effect = effect;
	}

	public void setExonId(String exonId) {
		this.exonId = exonId;
	}

	public void setFunClass(ChangeEffect.FunctionalClass funClass) {
		this.funClass = funClass;
	}

	public void setGene(String gene) {
		this.gene = gene;
	}

	public void setImpact(ChangeEffect.EffectImpact impact) {
		this.impact = impact;
	}

	public void setTranscriptId(String transcriptId) {
		this.transcriptId = transcriptId;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(effect);
		sb.append("(");

		if (impact != null) sb.append(impact);
		sb.append("|");

		if (funClass != null) sb.append(funClass);
		sb.append("|");

		sb.append(codon);
		sb.append("|");

		sb.append(aa);
		sb.append("|");

		if (aaLen > 0) sb.append(aaLen);
		sb.append("|");

		sb.append(gene);
		sb.append("|");

		sb.append(bioType);
		sb.append("|");

		if (coding != null) sb.append(coding);
		sb.append("|");

		sb.append(transcriptId);
		sb.append("|");

		sb.append(exonId);

		sb.append(")");

		return sb.toString();
	}
}
