package ca.mcgill.mcb.pcingola.interval;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;

import ca.mcgill.mcb.pcingola.interval.tree.IntervalForest;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * A collection of markers
 * 
 * @author pcingola
 *
 */
public class Markers implements Iterable<Marker>, Serializable {

	private static final long serialVersionUID = 259791388087691277L;
	boolean verbose = true;
	ArrayList<Marker> intervals;

	public Markers() {
		intervals = new ArrayList<Marker>();
	}

	/**
	 * Add an interval to the collection
	 * @param marker
	 */
	public void add(Marker marker) {
		intervals.add(marker);
	}

	/**
	 * Add all intervals
	 * @param intervalsMarkerIntervaloAdd
	 */
	public void add(Markers intervalsMarkerIntervaloAdd) {
		intervals.addAll(intervalsMarkerIntervaloAdd.intervals);
	}

	/**
	 * Are all intervals equal?
	 * @param intervals
	 * @return
	 */
	public boolean equals(Markers intervals) {
		if (intervals == null) return false;
		if (size() != intervals.size()) return false;

		// Sort both collections
		sort(false, false);
		intervals.sort(false, false);

		// Compare all intervals
		Iterator<Marker> it1 = iterator();
		Iterator<Marker> it2 = intervals.iterator();
		while (it1.hasNext() && it2.hasNext()) {
			Interval i1 = it1.next();
			Interval i2 = it2.next();
			if (!i1.equals(i2)) return false;
		}

		return true;
	}

	/**
	 * MarkerIntervalhis is a brute force implementation of 'intersect'.
	 * 
	 * WARNING: This method should only be used for debugging (or in very small collections) since it is extremely inefficient.
	 * For an efficient implementation, use IntervalForest class.
	 * 
	 * @param interval
	 * @return
	 */
	public Markers intersects(Marker interval) {
		Markers ints = new Markers();
		for (Marker i : this)
			if (i.intersects(interval)) ints.add(i);
		return ints;
	}

	public boolean isEmpty() {
		return intervals.isEmpty();
	}

	@Override
	public Iterator<Marker> iterator() {
		return intervals.iterator();
	}

	public Markers merge() {
		// Intervals sorted by start
		Markers intsSorted = new Markers();
		intsSorted.add(this);
		intsSorted.sort(false, false);

		// Merge intervals
		Markers intsMerged = new Markers();
		String tag = "", chromoName = "";
		Chromosome chromo = null;
		int start = -1, end = -1;
		for (Marker i : intsSorted) {

			// Different chromosome? => Re-start
			Chromosome ichromo = i.getChromosome();
			String ichromoName = ichromo.getId();
			if (!chromoName.equals(ichromoName)) {
				// Save current interval (if a any)
				if ((start >= 0) && (end >= 0)) {
					Marker im = new Marker(chromo, start, end, 1, tag);
					intsMerged.add(im);
				}

				chromoName = ichromoName;
				chromo = ichromo;
				start = end = -1;
				tag = "";
			}

			// Previous interval finished? => add it to list
			if (i.start > end) {
				if ((start >= 0) && (end >= 0)) {
					Marker im = new Marker(chromo, start, end, 1, tag);
					intsMerged.add(im);
				}
				start = end = -1;
				tag = "";
			}

			// Update interval 'start'
			if (start < 0) start = i.start;

			// Update 'end'
			end = Math.max(end, i.end);

			// Update tag
			if (tag.length() <= 0) tag = i.id;
			else tag += " " + i.id;
		}

		if ((start >= 0) && (end >= 0)) {
			Marker im = new Marker(chromo, start, end, 1, tag);
			intsMerged.add(im);
		}

		return intsMerged;
	}

	/**
	 * Calculate 'set minus' using one interval
	 * @param interval
	 * @return
	 */
	public Markers minus(Marker interval) {
		Markers ints = new Markers();

		// Add all intervals in 'this'
		for (Marker i : this)
			if (i.intersects(interval)) {
				if ((interval.getStart() <= i.getStart()) && (i.getEnd() <= interval.getEnd())) {
					// 'i' is included in 'interval' => Do not add 'i'
				} else if ((interval.getStart() <= i.getStart()) && (interval.getEnd() < i.getEnd())) {
					// 'interval' overlaps left part of 'i' => Include right part of 'i'
					ints.add(new Marker(i.getParent(), interval.getEnd() + 1, i.getEnd(), i.getStrand(), i.getId()));
				} else if ((i.getStart() < interval.getStart()) && (i.getEnd() <= interval.getEnd())) {
					// 'interval' overlaps right part of 'i' => Include left part of 'i'
					ints.add(new Marker(i.getParent(), i.getStart(), interval.getStart() - 1, i.getStrand(), i.getId()));
				} else if ((i.getStart() < interval.getStart()) && (interval.getEnd() < i.getEnd())) {
					// 'interval' overlaps middle of 'i' => Include left and right part of 'i'
					ints.add(new Marker(i.getParent(), i.getStart(), interval.getStart() - 1, i.getStrand(), i.getId()));
					ints.add(new Marker(i.getParent(), interval.getEnd() + 1, i.getEnd(), i.getStrand(), i.getId()));
				} else throw new RuntimeException("Interval intersection not analysed. This should nbever happen!");
			} else ints.add(i); // No intersection => Just add interval

		return ints;
	}

	/**
	 * Returns the result of this set minus 'intervals'
	 * 
	 * WARNING: This method should only be used for debugging (or in very small collections) since it is extremely inefficient.
	 * 
	 * @param interval
	 * @return
	 */
	public Markers minus(Markers intervals) {
		Markers result = new Markers();
		result.add(this);

		// Calculate 'set minus' for all 'intervals'
		for (Marker j : intervals)
			result = result.minus(j);

		return result;
	}

	/**
	 * Return a random interval within this collection
	 * @return
	 */
	public Interval rand() {
		int idx = (int) (Math.random() * intervals.size());
		return intervals.get(idx);
	}

	/**
	 * Read intervals from a file
	 * @param fileName
	 */
	public void read(String fileName, Genome genome) {
		// Read file
		String file = Gpr.readFile(fileName);

		// Parse lines
		String lines[] = file.split("\n");
		int lineNum = 1;
		for (String line : lines) {
			Marker interval = new Marker(null, 0, 0, 0, "");
			interval.parse(line, lineNum, genome, 0);
			add(interval);
			lineNum++;
		}
	}

	//	/**
	//	 * Save to a file
	//	 * @param fileName
	//	 */
	//	public void save(String fileName) {
	//		Gpr.toFile(fileName, toStringSave());
	//	}

	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}

	public int size() {
		return intervals.size();
	}

	/**
	 * Sort intervals
	 * @param byEnd : If true, sort by end. Otherwise sort by start
	 * @param reverse : Reverse order
	 */
	public void sort(boolean byEnd, boolean reverse) {
		if (byEnd) Collections.sort(intervals, new IntervalComparatorByEnd(reverse));
		else Collections.sort(intervals, new IntervalComparatorByStart(reverse));
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (Marker i : this)
			sb.append(i + " ");
		return sb.toString();
	}

	/**
	 * Show all intervals as an ASCII art
	 * @param maxLen
	 * @return
	 */
	public String toStringAsciiArt(int maxLen) {
		StringBuilder sb = new StringBuilder();

		// Separator
		String sep = "";
		for (int i = 0; i < maxLen; i++)
			sep += "=";

		// Show intervals
		String ch = "";
		for (Marker i : this) {
			if (!i.getChromosomeName().equals(ch)) {
				sb.append("|" + sep + "|\n");
				ch = i.getChromosomeName();
			}

			sb.append("|" + i.toStringAsciiArt(maxLen) + "|\t" + i.getChromosomeName() + ": [" + i.start + " - " + i.end + "] ");
			if ((i.id != null) && (i.id.length() > 0)) sb.append("'" + i.id + "'"); // Show tag (if any)
			sb.append("\n");
		}
		sb.append("|" + sep + "|\n");

		return sb.toString();
	}

	//	public String toStringSave() {
	//		StringBuilder sb = new StringBuilder();
	//		for (Marker i : this)
	//			sb.append(i.serializeSave() + "\n");
	//		return sb.toString();
	//	}

	/**
	 * Create a new 'intervals' object containing the union of both intervals
	 * @param intervalsById
	 */
	public Markers union(Markers intervalsForUnion) {
		Markers ints = new Markers();
		ints.add(this);
		ints.add(intervalsForUnion);
		return ints;
	}

	/**
	 * Perform the union of all overlapping intervals
	 * 
	 * For each marker, calculate all overlapping markers and create a new marker that contains them all.
	 * Return a set of those new markers.
	 * 
	 * @param markerIntervals
	 * @return
	 */
	public Markers unionOfOverlaps() {
		Markers unionOfOverlaps = new Markers();
		IntervalForest forest = new IntervalForest(this);

		HashSet<Marker> done = new HashSet<Marker>();
		for (Marker mi : this) {
			if (!done.contains(mi)) { // No added yet?
				Markers query = forest.query(mi);

				// Get union
				Marker union = new Marker(mi.getParent(), mi.getStart(), mi.getEnd(), mi.getStrand(), "");
				done.add(mi);
				for (Marker m : query) {
					if ((union.getStart() > m.getStart()) || (union.getEnd() < m.getEnd())) {
						// Create union result
						int start = Math.min(union.getStart(), m.getStart());
						int end = Math.max(union.getEnd(), m.getEnd());
						union = new Marker(mi.getParent(), start, end, m.getStrand(), "");
					}
					done.add(m);
				}

				// Add union
				unionOfOverlaps.add(union);
			}
		}

		return unionOfOverlaps;
	}
}
