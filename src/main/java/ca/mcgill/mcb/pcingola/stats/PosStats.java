package ca.mcgill.mcb.pcingola.stats;

import java.util.Random;

import ca.mcgill.mcb.pcingola.interval.Marker;

/**
 * How many changes per position do we have in a chromosome.
 * Summary by dividing the chromosome into MAX_BINS bins
 * 
 * @author pcingola
 */
public class PosStats extends ChrPosStats {

	public static final int DEFAULT_BINS = 100;

	int maxIndex = 0;

	public PosStats() {
		super("", DEFAULT_BINS);
		name = "";
		maxBins = DEFAULT_BINS;
		init(maxBins);
		factor = 1;
	}

	public PosStats(String name) {
		super(name, DEFAULT_BINS);
		this.name = name;
		maxBins = DEFAULT_BINS;
		init(maxBins);
		factor = 1;
	}

	public PosStats(String name, int maxBins) {
		super(name, maxBins);
		this.name = name;
		this.maxBins = maxBins;
		init(maxBins);
	}

	/**
	 * Create random counts (used for debugging)
	 * @param maxLen
	 * @param countMax
	 */
	public void rand(int maxLen, int countMax) {
		Random rand = new Random();
		for (int i = 0; i < Math.min(maxLen, count.length); i++)
			count[i] = rand.nextInt(countMax);
		maxIndex = maxLen;
	}

	/**
	 * Use 'num' as a sample
	 * @param num
	 */
	public void sample(Marker marker, Marker markerReference) {
		if (markerReference.intersects(marker)) {
			double step = 1;
			if (markerReference.size() > count.length) step = ((double) markerReference.size()) / count.length;

			int start = Math.max(marker.getStart(), markerReference.getStart());
			int end = Math.min(marker.getEnd(), markerReference.getEnd());

			double pos = start;
			int jmin = (int) ((start - markerReference.getStart()) / step);

			//			Gpr.debug("Name: "+ name + "Start: " + start + "\tend:" + end + "\tStep: " + step);

			int j;
			for (j = jmin; pos <= end; pos += step, j++)
				count[j]++;

			// Update maxIndex
			maxIndex = Math.max(maxIndex, j - 1);
		}
	}

	public int size() {
		return maxIndex;
	}

}
