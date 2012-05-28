package ca.mcgill.mcb.pcingola.nmer;

import gnu.trove.map.hash.TLongIntHashMap;
import gnu.trove.procedure.TLongIntProcedure;

import java.io.Serializable;

import ca.mcgill.mcb.pcingola.stats.Counter;

/**
 * Mark if an Nmer has been 'seen'
 * It only count up to 255 (one byte per counter)
 * 
 * @author pcingola
 */
public class NmerCount implements Serializable {

	private static final long serialVersionUID = 1L;
	public static boolean debug = false;

	int nmerSize;
	TLongIntHashMap hash;

	public NmerCount(int nmerSize) {
		this.nmerSize = nmerSize;
		hash = new TLongIntHashMap();
	}

	/**
	 * Average number of nmers
	 * @param threshold
	 * @return
	 */
	public double avg() {
		double total = total();
		double size = size();
		return size > 0 ? total / size : 0;
	}

	/**
	 * Count an instance of this Nmer
	 * @param nmer
	 */
	public void count(Nmer nmer) {
		long key = nmer.getNmer();
		int count = hash.get(key) + 1;
		hash.put(key, count);
	}

	/**
	 * Count how many nmers are below a given threshold
	 * @param threshold
	 * @return
	 */
	public long countLessThan(final int threshold) {
		final Counter counter = new Counter();
		hash.forEachEntry(new TLongIntProcedure() {

			@Override
			public boolean execute(long key, int value) {
				if (value > threshold) counter.inc();
				return true;
			}
		});
		return counter.count;
	}

	/**
	 * Max nmer count
	 * @param threshold
	 * @return
	 */
	public long max() {
		final Counter counter = new Counter();
		hash.forEachEntry(new TLongIntProcedure() {

			@Override
			public boolean execute(long key, int value) {
				if (value > counter.get()) counter.set(value);
				return true;
			}
		});
		return counter.count;
	}

	public int size() {
		return hash.size();
	}

	@Override
	public String toString() {
		return "Size: " + hash.size() + "\tTotal: " + total() + "\tAvg: " + avg() + "\tMax: " + max();
	}

	public String toStringAll() {
		return toStringAll(0);
	}

	public String toStringAll(final int minCount) {
		final StringBuilder sb = new StringBuilder();
		final Nmer nmer = new Nmer(nmerSize);
		hash.forEachEntry(new TLongIntProcedure() {

			@Override
			public boolean execute(long key, int value) {
				nmer.setNmer(key);
				if (value >= minCount) sb.append(nmer + "\t" + value + "\n");
				return true;
			}
		});
		return sb.toString();
	}

	/**
	 * Total number of nmers
	 * @param threshold
	 * @return
	 */
	public long total() {
		final Counter sum = new Counter();
		hash.forEachEntry(new TLongIntProcedure() {

			@Override
			public boolean execute(long key, int value) {
				sum.set(sum.get() + value);
				return true;
			}
		});
		return sum.get();
	}
}
