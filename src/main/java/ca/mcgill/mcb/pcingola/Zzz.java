package ca.mcgill.mcb.pcingola;

import ca.mcgill.mcb.pcingola.genBank.Feature;
import ca.mcgill.mcb.pcingola.genBank.GenBank;
import ca.mcgill.mcb.pcingola.util.Gpr;

public class Zzz {

	public static void main(String[] args) {
		GenBank gb = new GenBank(Gpr.HOME + "/genes.gb");
		for (Feature f : gb.getFeatures())
			System.out.println(f);
	}
}
