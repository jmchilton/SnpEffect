package ca.mcgill.mcb.pcingola.spliceSites;

import java.util.ArrayList;
import java.util.List;

/**
 * ACGT tree
 * 
 * @author pcingola
 */
public class AcgtTree {

	public static final char BASES[] = { 'A', 'C', 'G', 'T' };
	public static final double LOG2 = Math.log(2.0);
	public static final int FAKE_COUNTS = 1;
	public static final double MAX_ENTROPY = 2.0; // Maximum possible entropy for 4 symbols: - 4 * 1/4 * log2(1/4)

	String name;
	AcgtTree nodes[];
	int counts[];

	public static int base2index(char base) {
		switch (Character.toUpperCase(base)) {
		case 'A':
			return 0;
		case 'C':
			return 1;
		case 'G':
			return 2;
		case 'T':
			return 3;
		}
		throw new RuntimeException("Unknown base '" + base + "'");
	}

	public AcgtTree(String name) {
		this.name = name;
		nodes = new AcgtTree[4];
		counts = new int[4];
	}

	public void add(String sequence) {
		if ((sequence == null) || sequence.isEmpty()) return;

		// Count for this node
		char base = sequence.charAt(0);
		inc(base);
		AcgtTree node = getOrCreate(base);

		// Recurse into tree
		node.add(sequence.substring(1));
	}

	/**
	 * Calculate the entropy
	 * @return
	 */
	public double entropy() {
		double entropy = 0;
		for (double p : p())
			entropy += -p * Math.log(p) / LOG2;

		return entropy;
	}

	/**
	 * Find node names that are within the thresholds
	 * @param thresholdEntropy
	 * @param thresholdCount
	 * @return
	 */
	public List<String> findNodeNames(double thresholdEntropy, int thresholdCount) {
		ArrayList<String> names = new ArrayList<String>();
		if (totalCount() == 0) return names;

		names.add(name);

		for (char base : BASES) {
			int idx = base2index(base);
			AcgtTree n = nodes[idx];
			if ((n != null) && (n.entropy() <= thresholdEntropy) && (counts[idx] >= thresholdCount)) {
				names.addAll(n.findNodeNames(thresholdEntropy, thresholdCount));
			}
		}

		return names;
	}

	/**
	 * Get a node
	 * @param base
	 * @return
	 */
	public AcgtTree get(char base) {
		return nodes[base2index(base)];
	}

	/**
	 * Get node indexed by this string
	 * @param bases
	 * @return
	 */
	public AcgtTree get(String bases) {
		if (bases.isEmpty()) return this;
		char base = bases.charAt(0);
		AcgtTree node = get(base);
		if (node == null) return null;
		return node.get(bases.substring(1));
	}

	/**
	 * Get a node (create it if it doesn't exist)
	 * @param base
	 * @return
	 */
	public AcgtTree getOrCreate(char base) {
		AcgtTree node = get(base);
		if (node != null) return node;

		// Create node
		node = new AcgtTree(name + base);
		set(base, node);
		return node;
	}

	/**
	 * Increment counter for a base
	 * @param base
	 */
	public void inc(char base) {
		counts[base2index(base)]++;
	}

	double[] p() {
		int tot = 0;
		for (int c : counts) {
			c += FAKE_COUNTS;
			tot += c;
		}

		double p[] = new double[4];
		int i = 0;
		for (int c : counts) {
			c += FAKE_COUNTS;
			p[i] = ((double) c) / tot;
			i++;
		}

		return p;
	}

	public double seqConservation() {
		return (MAX_ENTROPY - entropy()) / MAX_ENTROPY;
	}

	/**
	 * Set a node
	 * @param base
	 * @param nodes
	 */
	public void set(char base, AcgtTree n) {
		nodes[base2index(base)] = n;
	}

	@Override
	public String toString() {
		return toString("", 2.0, 0);
	}

	public String toString(String tabs, double thresholdEntropy, int thresholdCount) {
		if (totalCount() == 0) return "";

		StringBuilder sb = new StringBuilder();
		for (char base : BASES) {
			int idx = base2index(base);
			AcgtTree n = nodes[idx];
			if (n != null) {
				sb.append(tabs + name + base + ": " + counts[idx] + "\t" + n.entropy() + "\n");
				if ((n.entropy() <= thresholdEntropy) && (counts[idx] >= thresholdCount)) sb.append(n.toString(tabs + "\t", thresholdEntropy, thresholdCount));
			}
		}

		return sb.toString();
	}

	public int totalCount() {
		int tot = 0;
		for (int c : counts)
			tot += c;
		return tot;
	}
}
