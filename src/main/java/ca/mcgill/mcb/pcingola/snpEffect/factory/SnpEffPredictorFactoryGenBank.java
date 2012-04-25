package ca.mcgill.mcb.pcingola.snpEffect.factory;

import ca.mcgill.mcb.pcingola.genBank.Features;
import ca.mcgill.mcb.pcingola.genBank.GenBank;
import ca.mcgill.mcb.pcingola.snpEffect.Config;

/**
 * This class creates a SnpEffectPredictor from a GenBank file.
 * 
 * @author pcingola
 */
public class SnpEffPredictorFactoryGenBank extends SnpEffPredictorFactoryFeatures {

	public SnpEffPredictorFactoryGenBank(Config config) {
		super(config);
		fileName = config.getBaseFileNameGenes() + ".gb";
	}

	/**
	 * Get features from file
	 * @return
	 */
	@Override
	protected Features readFeatures() {
		return new GenBank(fileName);
	}
}
