package ca.mcgill.mcb.pcingola.interval;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

/**
 * Interval that contains sub intervals.
 * 
 * @author pcingola
 *
 */
public class IntervalAndSubIntervals<T extends Marker> extends Marker implements Iterable<T> {

	private static final long serialVersionUID = 1636197649250882952L;
	HashMap<String, T> subIntervals;

	public IntervalAndSubIntervals() {
		super();
		subIntervals = new HashMap<String, T>();
	}

	public IntervalAndSubIntervals(Marker parent, int start, int end, int strand, String id) {
		super(parent, start, end, strand, id);
		subIntervals = new HashMap<String, T>();
	}

	/**
	 * Add a subinterval
	 * @param t
	 */
	public void add(T t) {
		if (subIntervals.put(t.getId(), t) != null) throw new RuntimeException(t.getClass().getSimpleName() + " '" + t.getId() + "' is already in " + this.getClass().getSimpleName() + " '" + id + "'");
	}

	/**
	 * Obtain a subinterval
	 * @param id
	 * @return
	 */
	public T get(String id) {
		return subIntervals.get(id);
	}

	@Override
	public Iterator<T> iterator() {
		return subIntervals.values().iterator();
	}

	public int numChilds() {
		return (subIntervals != null ? subIntervals.size() : 0);
	}

	/**
	 * Remove a subinterval
	 * @param t
	 */
	public void remove(T t) {
		subIntervals.remove(t.getId());
	}

	/**
	 * Parse a line from a serialized file
	 * @param line
	 * @return
	 */
	@SuppressWarnings("unchecked")
	@Override
	public void serializeParse(MarkerSerializer markerSerializer) {
		super.serializeParse(markerSerializer);

		Markers markers = markerSerializer.getNextFieldMarkers();
		for (Marker m : markers)
			add((T) m);
	}

	/**
	 * Create a string to serialize to a file
	 * @return
	 */
	@SuppressWarnings("unchecked")
	@Override
	public String serializeSave(MarkerSerializer markerSerializer) {
		return super.serializeSave(markerSerializer) + "\t" + markerSerializer.save((Collection<Marker>) subIntervals.values());
	}

	/**
	 * Return a collection of sub intervals sorted by natural order
	 * @return
	 */
	public List<T> sorted() {
		ArrayList<T> list = new ArrayList<T>();
		list.addAll(subIntervals.values());
		Collections.sort(list);
		return list;
	}

	/**
	 * Return a collection of sub intervals sorted by start position (if strand is >= 0) or 
	 * by reverse end position (if strans < 0) 
	 * @return
	 */
	public List<T> sortedStrand() {
		ArrayList<T> list = new ArrayList<T>();
		list.addAll(subIntervals.values());

		if (strand >= 0) Collections.sort(list, new IntervalComparatorByStart()); // Sort by start position 
		else Collections.sort(list, new IntervalComparatorByEnd(true)); // Sort by end position (reversed) 

		return list;
	}

	/**
	 * Return a collection of sub intervals
	 * @return
	 */
	public Collection<T> subintervals() {
		return subIntervals.values();
	}
}
