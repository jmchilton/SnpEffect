package ca.mcgill.mcb.pcingola;

import java.util.HashMap;
import java.util.Random;

import ca.mcgill.mcb.pcingola.stats.IntStats;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.GprSeq;

/**
 * Simple free energy calculation from stacked RNA constructs.
 * 
 * @author pablocingolani
 */
public class RnaStackedFreeEnergy {

	public static final int TO_INT = 1000;

	static String BASES[] = { "A", "C", "G", "U" };

	HashMap<String, Double> freeEnergy;

	public RnaStackedFreeEnergy(String detlaGfile) {
		freeEnergy = new HashMap<String, Double>();
		load(detlaGfile);
	}

	/**
	 * Calculate an empirical distribution
	 * @param size
	 * @param iterations
	 * @return
	 */
	public IntStats empiricalDistribution(int size, int iterations) {
		IntStats distribution = new IntStats();

		Random random = new Random();
		char bases1[] = new char[size];
		char bases2[] = new char[size];

		for (int i = 0; i < iterations; i++) {
			// Create random sequences
			for (int j = 0; j < bases1.length; j++) {
				bases1[j] = GprSeq.randBase(random);
				bases2[j] = GprSeq.randBase(random);
			}

			String seq1 = new String(bases1);
			String seq2 = new String(bases2);

			int energyInt = energyInt(seq1, seq2); // Get energy as an integer
			distribution.sample(energyInt); // Sample distrbiution
		}

		return distribution;
	}

	/** 
	 * Free energy (from table) 
	 * 		5' ==> 3'
	 * 		   WX
	 * 		   ZY
	 * 		3'<== 5'
	*/
	double energy(char w, char x, char y, char z) {
		String key = "" + w + x + y + z;
		Double d = freeEnergy.get(key);
		if (d == null) throw new RuntimeException("Unkonw WXYZ : '" + key + "'");
		return d;
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
	 * Return energy as an integer
	 * @param rna1
	 * @param rna2
	 * @return
	 */
	public int energyInt(String rna1, String rna2) {
		return (int) (energy(rna1, rna2) * TO_INT);
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