package ca.mcgill.mcb.pcingola;

import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.spliceSites.SpliceTypes;

public class Zzz {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Config config = new Config("amel2_cuff");
		//Config config = new Config("hg19");
		//		Config config = new Config("testHg3763Chr1");

		SpliceTypes spliceTypes = new SpliceTypes(config);
		spliceTypes.setVerbose(true);
		spliceTypes.analyzeAndCreate();
	}
}
