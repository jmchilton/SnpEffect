package ca.mcgill.mcb.pcingola.spliceSites;

/**
 * ACGT tree
 * 
 * @author pcingola
 */
public class AcgtTree {

	public static final char BASES[] = { 'A', 'C', 'G', 'T' };

	String name;
	AcgtTree node[];
	int count[];

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
		node = new AcgtTree[4];
		count = new int[4];
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
	 * Get a node
	 * @param base
	 * @return
	 */
	public AcgtTree get(char base) {
		return node[base2index(base)];
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
		count[base2index(base)]++;
	}

	/**
	 * Set a node
	 * @param base
	 * @param node
	 */
	public void set(char base, AcgtTree n) {
		node[base2index(base)] = n;
	}

	@Override
	public String toString() {
		return toString("", 0);
	}

	public String toString(String tabs, int threshold) {
		int total = 0;
		for (int c : count)
			total += c;
		if ((total) == 0) return "";

		StringBuilder sb = new StringBuilder();
		for (char base : BASES) {
			int idx = base2index(base);
			AcgtTree n = node[idx];
			sb.append(tabs + name + base + ": " + count[idx] + "\n");
			if ((n != null) && (count[idx] >= threshold)) sb.append(n.toString(tabs + "\t", threshold));
		}

		return sb.toString();
	}

}
