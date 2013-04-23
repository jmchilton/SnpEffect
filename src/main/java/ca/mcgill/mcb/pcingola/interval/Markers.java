package ca.mcgill.mcb.pcingola.interval;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;

import ca.mcgill.mcb.pcingola.fileIterator.BedFileIterator;
import ca.mcgill.mcb.pcingola.fileIterator.BigBedFileIterator;
import ca.mcgill.mcb.pcingola.interval.tree.IntervalForest;
import ca.mcgill.mcb.pcingola.serializer.MarkerSerializer;

/**
 * A collection of markers
 * 
 * @author pcingola
 */
public class Markers implements Serializable, Collection<Marker> {

	private static final long serialVersionUID = 259791388087691277L;
	protected ArrayList<Marker> markers;
	protected String name = "";

	/**
	 * Read markers from a file
	 * Supported formats: BED, TXT, BigBed
	 */
	public static Markers readMarkers(String fileName) {
		String flLower = fileName.toLowerCase();

		// Load according to file type
		if (flLower.endsWith(".txt") || flLower.endsWith(".txt.gz")) return new BedFileIterator(fileName, null, 1).loadMarkers(); // TXT is assumed to be "chr \t start \t end"
		else if (flLower.endsWith(".bed") || flLower.endsWith(".bed.gz")) return new BedFileIterator(fileName).loadMarkers();
		else if (flLower.endsWith(".bb")) return new BigBedFileIterator(fileName).loadMarkers();
		else throw new RuntimeException("Unrecognized genomig interval file type '" + fileName + "'");
	}

	public Markers() {
		markers = new ArrayList<Marker>();
	}

	@SuppressWarnings({ "rawtypes", "unchecked" })
	public Markers(Collection otherMarkers) {
		markers = new ArrayList<Marker>();
		addAll(otherMarkers);
	}

	public Markers(Markers otherMarkers) {
		markers = new ArrayList<Marker>();
		addAll(otherMarkers.getMarkers());
	}

	public Markers(String name) {
		this.name = name;
		markers = new ArrayList<Marker>();
	}

	/**
	 * Add an interval to the collection
	 * @param marker
	 */
	@Override
	public boolean add(Marker marker) {
		return markers.add(marker);
	}

	/**
	 * Add all intervals
	 * @param markersToAdd
	 */
	public Markers add(Markers markersToAdd) {
		markers.addAll(markersToAdd.markers);
		return this;
	}

	/**
	 * Add all markers in this collection
	 * @param intervalsMarkerIntervaloAdd
	 */
	@Override
	public boolean addAll(Collection<? extends Marker> mm) {
		boolean changed = false;
		for (Marker m : mm)
			changed |= markers.add(m);
		return changed;
	}

	@Override
	public void clear() {
		markers.clear();
	}

	@Override
	public boolean contains(Object o) {
		return markers.contains(o);
	}

	@Override
	public boolean containsAll(Collection<?> c) {
		return markers.containsAll(c);
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

	public List<Marker> getMarkers() {
		return markers;
	}

	public String getName() {
		return name;
	}

	/**
	 * Perform the intersection of all overlapping intervals
	 * 
	 * For each marker, calculate all overlapping markers and create a new marker that contains them all.
	 * Return a set of those new markers.
	 * 
	 * @param markerIntervals
	 * @return
	 */
	public Markers intersect() {
		Markers intersectOfOverlaps = new Markers();
		IntervalForest forest = new IntervalForest(this);

		HashSet<Marker> done = new HashSet<Marker>();
		for (Marker mi : this) {
			if (!done.contains(mi)) { // No added yet?
				Markers query = forest.query(mi);

				// Get intersect
				Marker intersect = new Marker(mi.getParent(), mi.getStart(), mi.getEnd(), mi.getStrand(), "");
				done.add(mi);
				for (Marker m : query) {
					if (intersect != null) {
						if ((intersect.getStart() < m.getStart()) || (intersect.getEnd() > m.getEnd())) {
							intersect = intersect.intersect(m);
						}
					}
					done.add(m);
				}

				// Add union
				if (intersect != null) intersectOfOverlaps.add(intersect);
			}
		}

		return intersectOfOverlaps;
	}

	@Override
	public boolean isEmpty() {
		return markers.isEmpty();
	}

	@Override
	public Iterator<Marker> iterator() {
		return markers.iterator();
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
		int idx = (int) (Math.random() * markers.size());
		return markers.get(idx);
	}

	@Override
	public boolean remove(Object o) {
		return markers.remove(o);
	}

	@Override
	public boolean removeAll(Collection<?> c) {
		return markers.removeAll(c);
	}

	@Override
	public boolean retainAll(Collection<?> c) {
		return markers.retainAll(c);
	}

	/**
	 * Save to a file using a serializer
	 * @param fileName
	 */
	public void save(String fileName) {
		// Nothing to save
		if (size() <= 0) return;

		// We must add genome and all chromosomes to the list (otherwise the serializer cannot save all references)
		Genome genome = markers.get(0).getGenome();

		// Add all chromosomes to a set (to make sure they are added only once)
		HashSet<Chromosome> chromos = new HashSet<Chromosome>();
		for (Marker m : this)
			chromos.add(m.getChromosome());

		// Create a set of all markers to be saved
		Markers markersToSave = new Markers();

		// Add genome
		markersToSave.add(genome);

		// Add chromosomes
		for (Chromosome chr : chromos)
			markersToSave.add(chr);

		// Add markers
		for (Marker m : markers)
			markersToSave.add(m);

		// Save
		MarkerSerializer markerSerializer = new MarkerSerializer();
		markerSerializer.save(fileName, markersToSave);
	}

	@Override
	public int size() {
		return markers.size();
	}

	/**
	 * Sort intervals
	 */
	public Markers sort() {
		return sort(false, false);
	}

	/**
	 * Sort intervals
	 * @param byEnd : If true, sort by end. Otherwise sort by start
	 * @param reverse : Reverse order
	 */
	public Markers sort(boolean byEnd, boolean reverse) {
		if (byEnd) Collections.sort(markers, new IntervalComparatorByEnd(reverse));
		else Collections.sort(markers, new IntervalComparatorByStart(reverse));
		return this;
	}

	@Override
	public Object[] toArray() {
		return markers.toArray();
	}

	@Override
	public <T> T[] toArray(T[] a) {
		return markers.toArray(a);
	}

	@Override
	public String toString() {
		int num = 1;
		StringBuilder sb = new StringBuilder();
		for (Marker i : this)
			sb.append("\t" + (num++) + ":\t" + i.getChromosomeName() + "\t" + i.getStart() + "\t" + i.getEnd() + "\t" + i.getClass().getSimpleName() + "\t" + i.getId() + "\n");
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

	public String toStringTxt() {
		StringBuilder sb = new StringBuilder();
		for (Marker i : this)
			sb.append(i.getChromosomeName() + "\t" + i.getStart() + "\t" + i.getEnd() + "\t" + i.getId() + "\n");
		return sb.toString();
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
	public Markers union() {
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
					if ((union != null) && (union.getStart() > m.getStart()) || (union.getEnd() < m.getEnd())) union = union.union(m);
					done.add(m);
				}

				// Add union
				if (union != null) unionOfOverlaps.add(union);
			}
		}

		return unionOfOverlaps;
	}

	/**
	 * Remove duplicated markers
	 * @return this object 
	 */
	public Markers unique() {
		HashSet<Marker> set = new HashSet<Marker>();
		set.addAll(markers);
		markers = new ArrayList<Marker>();
		markers.addAll(set);
		return this;
	}

}
