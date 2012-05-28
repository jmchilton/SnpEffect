package ca.mcgill.mcb.pcingola;

import java.util.HashMap;

import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * Simple free energy calculation from stacked RNA constructs.
 * 
 * @author pablocingolani
 */
public class RnaStackedFreeEnergy {

	static String BASES[] = { "A", "C", "G", "U" };

	HashMap<String, Double> freeEnergy;

	public RnaStackedFreeEnergy(String detlaGfile) {
		freeEnergy = new HashMap<String, Double>();
		load(detlaGfile);
	}

	/** 
	 * Free energy (from table) 
	 * 		5' ==> 3'
	 * 		   WX
	 * 		   ZY
	 * 		3'<== 5'
	*/
	double energy(char w, char x, char y, char z) {
		return freeEnergy.get("" + w + x + y + z);
	}

	/**
	 * Calculate a simple approximation of free energy of two RNA strings
	 * 
	 * @param rna1
	 * @param rna2
	 * @return
	 */
	public double energy(String rna1, String rna2) {
		if (rna1.length() != rna2.length()) throw new RuntimeException("Only RNA with same length are supported now.");

		// Convert to upper case, convert 'T' to 'U'
		rna1 = rna1.toUpperCase().replace('T', 'U');
		rna2 = rna2.toUpperCase().replace('T', 'U');
		char rnaBases1[] = rna1.toCharArray();
		char rnaBases2[] = rna2.toCharArray();

		// Calculate free energy approximation
		int max = rnaBases1.length - 1;
		double energy = 0;
		for (int i = 0; i < max; i++)
			energy += energy(rnaBases1[i], rnaBases1[i + 1], rnaBases2[i + 1], rnaBases2[i]);

		return energy;
	}

	/**
	 * Load free energy values from a file
	 * @param detlaGfile
	 */
	void load(String detlaGfile) {
		String file = Gpr.readFile(detlaGfile);
		String lines[] = file.split("\n");

		// Parse each line
		for (String line : lines) {
			String recs[] = line.split("\\s");

			// Bases
			int rec = 0;
			String w = recs[rec++];
			String z = recs[rec++];

			// Matrix in one line
			for (String x : BASES) {
				for (String y : BASES) {
					String key = w + x + y + z;
					if (!recs[rec].equals(".")) freeEnergy.put(key, Gpr.parseDoubleSafe(recs[rec]));
					else freeEnergy.put(key, 0.0);
					rec++;
				}
			}
		}
	}
}