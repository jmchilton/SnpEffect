package ca.mcgill.mcb.pcingola;

import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;
import ca.mcgill.mcb.pcingola.vcf.VcfGenotype;

/**
 * Simple test program
 * @author pcingola
 */
public class ZzzVcfLoad {

	String vcfFileName;
	GenotypeVector genotypeVectors[];

	public static void main(String[] args) {
		ZzzVcfLoad vcfLoad = new ZzzVcfLoad();
		vcfLoad.parse(args);
		vcfLoad.run();
	}

	public ZzzVcfLoad() {
	}

	/**
	 * Parse command line arguments
	 * @param args
	 */
	public void parse(String[] args) {
		if (args.length != 1) {
			System.err.println("Usage: " + ZzzVcfLoad.class.getSimpleName() + " vcfFile");
			System.exit(-1);
		}

		vcfFileName = args[0];
	}

	public boolean run() {
		//---
		// Create data structure
		//---
		Timer.showStdErr("Counting lines form file" + vcfFileName);
		int numLines = Gpr.countLines(vcfFileName);
		Timer.showStdErr("Done. Number of lines: " + numLines);

		Timer.showStdErr("Loading file " + vcfFileName);
		VcfFileIterator vcf = new VcfFileIterator(vcfFileName);
		for (VcfEntry ve : vcf) {
			if (genotypeVectors == null) {
				Timer.showStdErr("Initializing data structures.");

				genotypeVectors = new GenotypeVector[ve.getVcfGenotypes().size()];
				for (int i = 0; i < genotypeVectors.length; i++)
					genotypeVectors[i] = new GenotypeVector(numLines);

				Timer.showStdErr("Done.");
			}

			//System.out.print(ve.getChromosomeName() + ":" + ve.getStart());
			for (VcfGenotype vg : ve) {
				int code = vg.getGenotypeCode();
				if (code < 0) code = 0;
				//System.out.print(String.format("%2d", code));
			}
			//System.out.println("");
		}

		Timer.showStdErr("Done");
		return true;
	}
}
