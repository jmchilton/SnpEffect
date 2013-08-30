package ca.mcgill.mcb.pcingola.reactome.events;

import java.util.HashSet;
import java.util.Iterator;

import ca.mcgill.mcb.pcingola.reactome.Entity;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * A Reactome complex (a bunch of molecules or complexes
 * 
 * @author pcingola
 *
 */
public class Complex extends Reaction implements Iterable<Entity> {

	public Complex(int id, String name) {
		super(id, name);
	}

	public void add(Entity e) {
		addInput(e);
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

		if (!Double.isNaN(fixedOutput)) output = fixedOutput;
		else {
			if (doneEntities.contains(this)) return output; // Make sure we don't calculate twice
			doneEntities.add(this); // Keep 'entities' set up to date

			// Sum inputs
			double sum = 0;
			int count = 0;
			for (Entity e : this) {
				double out = e.calc(doneEntities);
				if (e.hasOutput()) {
					count += inputs.get(e);
					sum += out * out;
				}
			}

			// Calculate output
			if (count > 0) output = cap(sum / Math.sqrt(count));
			else output = Double.NaN;
		}

		if (debug) System.out.println(output + "\tfixed:" + isFixed() + "\tid:" + id + "\ttype:" + getClass().getSimpleName() + "\tname:" + name);
		return output;
	}

	public boolean contains(Object o) {
		return inputs.containsKey(o);
	}

	public boolean isEmpty() {
		return inputs.isEmpty();
	}

	@Override
	public Iterator<Entity> iterator() {
		return inputs.keySet().iterator();
	}

	public int size() {
		return inputs.size();
	}

	@Override
	public String toString(int tabs, HashSet<Entity> done) {
		done.add(this);

		StringBuilder sb = new StringBuilder();
		sb.append(Gpr.tabs(tabs) + toStringSimple() + "\n");

		for (Entity e : this)
			sb.append(e.toString(tabs + 1, done) + "\n");

		return sb.toString();
	}

}
