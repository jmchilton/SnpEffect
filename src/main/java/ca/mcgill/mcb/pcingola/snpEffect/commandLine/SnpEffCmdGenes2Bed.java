package ca.mcgill.mcb.pcingola.snpEffect.commandLine;

import java.util.HashSet;

import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Genome;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * Simple test program
 * @author pcingola
 */
public class SnpEffCmdGenes2Bed extends SnpEff {

	HashSet<String> geneIds;

	public static void main(String[] args) {
		SnpEffCmdGenes2Bed conf2down = new SnpEffCmdGenes2Bed();
		conf2down.parseArgs(args);
		conf2down.run();
	}

	public SnpEffCmdGenes2Bed() {
		geneIds = new HashSet<String>();
	}

	@Override
	public void parseArgs(String[] args) {
		if (args.length < 1) usage(null);

		// Parse command line arguments
		genomeVer = args[0];
		for (int i = 1; i < args.length; i++)
			geneIds.add(args[i]);
	}

	@Override
	public boolean run() {
		// Load config & database
		if (verbose) Timer.showStdErr("Loading config file '" + configFile + "'");
		config = new Config(genomeVer, configFile);
		if (verbose) Timer.showStdErr("Loading database " + genomeVer);
		config.loadSnpEffectPredictor();
		Genome genome = config.getGenome();

		// Find genes
		if (verbose) Timer.showStdErr("Calculate");
		System.out.println("#chr\tstart\tend\tgeneName;geneId");
		for (Gene g : genome.getGenes()) {
			// Is gene.id or gene.name in geneSet? => Show it
			if (geneIds.contains(g.getId()) || geneIds.contains(g.getGeneName())) System.out.println(g.getChromosomeName() + "\t" + g.getStart() + "\t" + g.getEnd() + "\t" + g.getGeneName() + ";" + g.getId());
		}

		if (verbose) Timer.showStdErr("Done");

		return true;
	}

	@Override
	public void usage(String message) {
		if (message != null) System.err.println("Error: " + message + "\n");
		System.err.println("Usage: " + SnpEffCmdGenes2Bed.class.getSimpleName() + " genomeVer geneId_1 geneId_2 geneId_3 ... geneId_N");
		System.exit(-1);
	}

}
