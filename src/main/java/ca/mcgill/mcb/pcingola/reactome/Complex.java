package ca.mcgill.mcb.pcingola.reactome;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * A Reactome complex (a bunch of molecules or complexes
 * 
 * @author pcingola
 *
 */
public class Complex extends Entity implements Iterable<Entity> {

	HashMap<Entity, Integer> entities; // Count number of entities

	public Complex(int id, String name) {
		super(id, name);
		entities = new HashMap<Entity, Integer>();
	}

	public void add(Entity e) {
		if (e == null) return;

		if (!entities.containsKey(e)) entities.put(e, 1);
		entities.put(e, entities.get(e) + 1);
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

		doneEntities.add(this); // Keep 'entities' set up to date

		for (Entity e : this) {
			int count = entities.get(e);
			double out = e.calc(doneEntities);
			addWeight(out * count);
		}

		output = getWeight(); // Calculate output
		if (hasOutput()) Gpr.debug("COMPLEX:\t" + this.getName() + "\t=>\t" + output);

		return output;
	}

	public boolean contains(Object o) {
		return entities.containsKey(o);
	}

	public boolean isEmpty() {
		return entities.isEmpty();
	}

	@Override
	public Iterator<Entity> iterator() {
		return entities.keySet().iterator();
	}

	public int size() {
		return entities.size();
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(getClass().getSimpleName() + "[" + id + "]: " + name);

		for (Entity e : this)
			sb.append(e + " ");
		return sb.toString();
	}

}
