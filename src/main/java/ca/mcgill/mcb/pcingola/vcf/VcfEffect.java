package ca.mcgill.mcb.pcingola.vcf;

import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * An 'EFF' entry in a vcf line
 * 
 * @author pablocingolani
 */
public class VcfEffect {

	public static final boolean USE_OLD_FORMAT = true;

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

		try {
			// Parse each sub field
			int index = 0;
			if (effs[index].startsWith("REGULATION")) effect = ChangeEffect.EffectType.REGULATION;
			else effect = ChangeEffect.EffectType.valueOf(effs[index]);
			index++;

			if ((effs.length > index) && !effs[index].isEmpty()) impact = ChangeEffect.EffectImpact.valueOf(effs[index]);
			index++;

			if ((effs.length > index) && !effs[index].isEmpty()) funClass = ChangeEffect.FunctionalClass.valueOf(effs[index]);
			index++;

			if ((effs.length > index) && !effs[index].isEmpty()) codon = effs[index];
			index++;

			if ((effs.length > index) && !effs[index].isEmpty()) aa = effs[index];
			index++;

			if (!USE_OLD_FORMAT) {
				if ((effs.length > index) && !effs[index].isEmpty()) aaLen = Gpr.parseIntSafe(effs[index]);
				else aaLen = 0;
				index++;
			}

			if ((effs.length > index) && !effs[index].isEmpty()) gene = effs[index];
			index++;

			if ((effs.length > index) && !effs[index].isEmpty()) bioType = effs[index];
			index++;

			if ((effs.length > index) && !effs[index].isEmpty()) coding = ChangeEffect.Coding.valueOf(effs[index]);
			index++;

			if ((effs.length > index) && !effs[index].isEmpty()) transcriptId = effs[index];
			index++;

			if ((effs.length > index) && !effs[index].isEmpty()) exonId = effs[index];
			index++;
		} catch (Exception e) {
			String fields = "";
			for (int i = 0; i < effs.length; i++)
				fields += "\t" + i + " : '" + effs[i] + "'\n";
			throw new RuntimeException("Error parsing: '" + eff + "'\n" + fields, e);
		}
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
