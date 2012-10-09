package ca.mcgill.mcb.pcingola.snpEffect.factory;

import java.util.ArrayList;
import java.util.List;

import ca.mcgill.mcb.pcingola.genBank.Embl;
import ca.mcgill.mcb.pcingola.genBank.Features;
import ca.mcgill.mcb.pcingola.genBank.GenBank;
import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * This class creates a SnpEffectPredictor from an Embl file.
 * 
 * @author pcingola
 */
public class SnpEffPredictorFactoryEmbl extends SnpEffPredictorFactoryFeatures {

	public SnpEffPredictorFactoryEmbl(Config config) {
		super(config);
		fileName = config.getBaseFileNameGenes() + ".embl";
	}

	/**
	 * Get features from file
	 * @return
	 */
	@Override
	protected List<Features> readFeatures() {
		ArrayList<Features> featList = new ArrayList<Features>();

		if (Gpr.canRead(fileName)) {
			System.out.println("Reading data file  : '" + fileName + "'");
			featList.add(new Embl(fileName));
		} else {
			for (Chromosome chr : config.getGenome()) {
				String chrFileName = config.getDirDataVersion() + "/" + chr.getId() + ".embl";
				System.out.println("Reading data file  : '" + chrFileName + "'");
				featList.add(new GenBank(chrFileName));
			}
		}

		return featList;
	}
}
