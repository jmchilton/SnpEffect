package ca.mcgill.mcb.pcingola.reactome.events;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import ca.mcgill.mcb.pcingola.reactome.Entity;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * A Reactome pathway
 * 
 * @author pcingola
 *
 */
public class Pathway extends Reaction implements Iterable<Event> {

	HashMap<Event, Integer> events; // Count number of events

	public Pathway(int id, String name) {
		super(id, name);
		events = new HashMap<Event, Integer>();
	}

	public void add(Event e) {
		if (e == null) return;

		if (!events.containsKey(e)) events.put(e, 1);
		events.put(e, events.get(e) + 1);
	}

	/**
	 * Calculate entities.
	 * Make sure we don't calculate twice (keep 'doneEntities' set up to date)
	 * 
	 * @param doneEntities
	 * @return
	 */
	@Override
	public double calc(HashSet<Entity> doneEntities) {
		if (doneEntities.contains(this)) return output; // Make sure we don't calculate twice
		if (hasOutput()) return output; // TODO: Remove this!

		doneEntities.add(this); // Keep 'entities' set up to date

		for (Entity e : this) {
			double out = e.calc(doneEntities);
			addWeight(out);
		}

		output = getWeight();
		if (hasOutput()) Gpr.debug("PATHWAY:\t" + this.getName() + "\t=>\t" + output);

		return output;
	}

	public boolean contains(Object o) {
		return events.containsKey(o);
	}

	public boolean isEmpty() {
		return events.isEmpty();
	}

	public Iterator<Event> iterator() {
		return events.keySet().iterator();
	}

	public int size() {
		return events.size();
	}

	@Override
	public String toString() {
		return toString(0, new HashSet<Entity>());
	}

	@Override
	public String toString(int tabs, HashSet<Entity> done) {
		StringBuilder sb = new StringBuilder();
		sb.append(super.toString(tabs, done));

		for (Event e : this)
			sb.append(e.toString(tabs + 1, done) + "\n");

		return sb.toString();
	}

}
