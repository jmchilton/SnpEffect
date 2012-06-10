package ca.mcgill.mcb.pcingola.motif;

/**
 * Create a DNA motif count matrix
 * 
 * @author pcingola
 */
public class Pwm {

	public static final char BASES[] = { 'A', 'C', 'G', 'T' };
	int countMatrix[][]; // Keep counts for each base and position: countMatrix[base][position]
	int length;
	int totalCount;

	public Pwm(int length) {
		this.length = length;
		countMatrix = new int[4][length];
	}

	/**
	 * Transform a base into a code
	 * @param base
	 * @return
	 */
	int base2int(char base) {
		switch (base) {
		case 'a':
		case 'A':
			return 0;
		case 'c':
		case 'C':
			return 1;
		case 'g':
		case 'G':
			return 2;
		case 't':
		case 'T':
		case 'u':
		case 'U':
			return 3;
		}

		return -1;
	}

	/**
	 * Get counts for a given position
	 * @param base
	 * @param position
	 * @return
	 */
	public int getCount(char base, int position) {
		return countMatrix[base2int(base)][position];
	}

	public int getTotalCount() {
		return totalCount;
	}

	/**
	 * Matrix size
	 * @return
	 */
	int size() {
		return countMatrix[0].length;
	}

	@Override
	public String toString() {
		StringBuffer sb = new StringBuffer();
		for (int b = 0; b < BASES.length; b++) {
			sb.append(BASES[b] + "\t");
			for (int i = 0; i < countMatrix[b].length; i++)
				sb.append(countMatrix[b][i] + "\t");
			sb.append("\n");
		}

		sb.append("Max:\t");
		for (int i = 0; i < countMatrix[0].length; i++) {
			int max = 0, maxb = 0;
			for (int b = 0; b < BASES.length; b++) {
				if (max < countMatrix[b][i]) {
					max = countMatrix[b][i];
					maxb = b;
				}
			}
			sb.append(BASES[maxb] + "\t");
		}
		sb.append("\n");

		return sb.toString();
	}

	public void updateCounts(String dna) {
		updateCounts(dna, 1);
	}

	/**
	 * Update counts matrix.
	 * @param dna
	 */
	public void updateCounts(String dna, int inc) {
		totalCount += inc;
		char bases[] = dna.toCharArray();

		for (int i = 0; i < bases.length; i++) {
			int code = base2int(bases[i]);
			if (code >= 0) countMatrix[code][i] += inc;
		}
	}
}
