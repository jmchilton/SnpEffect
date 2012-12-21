package ca.mcgill.mcb.pcingola;

/**
 * A vector of genotypes in a 'compact' structure
 * 
 * Note: Genotypes 0/0, 0/1, 1/0, 1/1 are stored in 2 bits.
 * WARNIGN: Other genotypes are ignored or silently converted to 0/0
 *  
 * @author pcingola
 */
public class GenotypeVector {

	int size; // Size in elements (genotypes)
	byte genotype[];

	public GenotypeVector(int size) {
		this.size = size;
		genotype = new byte[pos2byte(size) + 1];

		for (int i = 0; i < genotype.length; i++)
			genotype[i] = 0;
	}

	/**
	 * Convert position to byte number
	 * @param pos
	 * @return
	 */
	int pos2byte(int pos) {
		return pos / 4 + ((pos & 0x03) == 0 ? 0 : 1);
	}

	public int size() {
		return size;
	}
}
