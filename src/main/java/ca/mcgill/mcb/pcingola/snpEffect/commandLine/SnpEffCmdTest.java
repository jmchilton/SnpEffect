package ca.mcgill.mcb.pcingola.snpEffect.commandLine;

import ca.mcgill.mcb.pcingola.gtex.Gtex;
import ca.mcgill.mcb.pcingola.gtex.GtexExperiment;
import ca.mcgill.mcb.pcingola.reactome.Reactome;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * Command line: Test
 * 
 * Note: Used for testing weird stuff
 * 
 * @author pcingola
 */
public class SnpEffCmdTest extends SnpEff {

	String reactomeDir = Gpr.HOME + "/snpEff/db/reactome/txt/";
	String geneIdsFile = Gpr.HOME + "/snpEff/db/reactome/gene_ids/biomart_query_uniq.txt";

	public SnpEffCmdTest() {
		super();
	}

	/**
	 * Parse command line arguments
	 */
	@Override
	public void parseArgs(String[] args) {
		//		if (args.length <4) usage(null);
	}

	/**
	 * Run command
	 */
	@Override
	public boolean run() {

		Reactome reactome = new Reactome(reactomeDir);
		reactome.load();
		reactome.loadGeneIds(geneIdsFile); // Load Gene IDs data

		//---
		// Load GTEX data
		//---
		Timer.showStdErr("Loading GTEx data");
		String gtexDir = Gpr.HOME + "/snpEff/db/GTEx";
		String gtexSamples = gtexDir + "/GTEx_Analysis_Annotations_Sample_DS__Pilot_2013_01_31.txt";
		String gtexData = gtexDir + "/gtex_norm.zzz.txt";
		Gtex gtex = new Gtex(gtexSamples, gtexData);
		gtex.getClass();

		//---
		// Simulate...!?
		//---
		for (GtexExperiment gtexExperiment : gtex) {
			if (gtexExperiment.size() > 0) reactome.zzz(gtexExperiment);
		}

		return true;
	}

	/**
	 * Show usage and exit
	 */
	@Override
	public void usage(String message) {
		if (message != null) System.err.println("Error: " + message + "\n");
		System.err.println("snpEff version " + SnpEff.VERSION);
		System.err.println("Usage: snpEff test [options] ...");
		System.exit(-1);
	}

}
