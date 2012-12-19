package ca.mcgill.mcb.pcingola.snpEffect.commandLine;

import java.util.HashSet;

import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Genome;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * Simple test program
 * @author pcingola
 */
public class SnpEffCmdGenes2Bed extends SnpEff {

	HashSet<String> geneIds;
	String fileName = null;

	public static void main(String[] args) {
		SnpEffCmdGenes2Bed conf2down = new SnpEffCmdGenes2Bed();
		conf2down.parseArgs(args);
		conf2down.run();
	}

	public SnpEffCmdGenes2Bed() {
		super();
		geneIds = new HashSet<String>();
	}

	void load() {
		// Parse file
		if (fileName != null) {
			if (verbose) Timer.showStdErr("Loading genes list from file '" + fileName + "'");

			String lines[] = Gpr.readFile(fileName).split("\n");
			if (lines.length <= 0) throw new RuntimeException("Cannot read file '" + fileName + "'");

			for (String line : lines) {
				String id = line.trim();
				if (!id.isEmpty()) geneIds.add(id);
			}
		}
	}

	@Override
	public void parseArgs(String[] args) {
		if (args.length < 1) usage(null);

		// Parse command line arguments
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-f")) { // List in a file?
				if ((i + 1) < args.length) fileName = args[++i];
				else usage("Option '-f' without file argument");
			} else if ((genomeVer == null) || genomeVer.isEmpty()) genomeVer = args[i];
			else geneIds.add(args[i]);
		}
	}

	@Override
	public boolean run() {
		load();
		if (verbose) Timer.showStdErr("Number of gene IDs to look up: " + geneIds.size());

		// Load config & database
		if (verbose) Timer.showStdErr("Loading config file '" + configFile + "'");
		config = new Config(genomeVer, configFile);
		if (verbose) Timer.showStdErr("Loading database " + genomeVer);
		config.loadSnpEffectPredictor();
		Genome genome = config.getGenome();

		// Find genes
		if (verbose) Timer.showStdErr("Finding genes.");
		int found = 0;
		System.out.println("#chr\tstart\tend\tgeneName;geneId");
		for (Gene g : genome.getGenes()) {
			// Is gene.id or gene.name in geneSet? => Show it
			if (geneIds.contains(g.getId()) || geneIds.contains(g.getGeneName())) {
				found++;
				System.out.println(g.getChromosomeName() + "\t" + g.getStart() + "\t" + g.getEnd() + "\t" + g.getGeneName() + ";" + g.getId());
			}
		}

		if (verbose) Timer.showStdErr("Done. Found " + found + " / " + geneIds.size());

		return true;
	}

	@Override
	public void usage(String message) {
		if (message != null) System.err.println("Error: " + message + "\n");
		System.err.println("Usage: " + SnpEffCmdGenes2Bed.class.getSimpleName() + " genomeVer [-f genes.txt | geneList]}");
		System.err.println("Options: ");
		System.err.println("\t-f <file.txt>  : A TXT file having one gene ID (or name) per line.");
		System.err.println("\tgeneList       : A list of gene IDs or names. One per command line argument: geneId_1 geneId_2 geneId_3 ... geneId_N");
		System.exit(-1);
	}

}
