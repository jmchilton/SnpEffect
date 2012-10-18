package ca.mcgill.mcb.pcingola;

import ca.mcgill.mcb.pcingola.genBank.EmblFile;
import ca.mcgill.mcb.pcingola.genBank.Features;
import ca.mcgill.mcb.pcingola.util.Gpr;

public class Zzz {

	public static void main(String[] args) {
		String file = Gpr.HOME + "/snpEff/data/spombe/genes.embl";

		EmblFile emblFile = new EmblFile(file);
		for (Features f : emblFile) {
			System.out.println("Reaging!");
		}

	}
}
