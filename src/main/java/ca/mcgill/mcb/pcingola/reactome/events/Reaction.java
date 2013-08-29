package ca.mcgill.mcb.pcingola.reactome.events;

import java.util.HashMap;
import java.util.HashSet;

import ca.mcgill.mcb.pcingola.reactome.Entity;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * A reaction
 * 
 * @author pcingola
 *
 */
public class Reaction extends Event {

	/**
	 * Reaction regulation types
	 * @author pablocingolani
	 *
	 */
	public enum RegulationType {
		NegativeRegulation, PositiveRegulation, Requirement
	};

	HashSet<Entity> inputs, outputs, catalyst;
	HashMap<Entity, RegulationType> regulator;

	public Reaction(int id, String name) {
		super(id, name);
		inputs = new HashSet<Entity>();
		outputs = new HashSet<Entity>();
		catalyst = new HashSet<Entity>();
		regulator = new HashMap<Entity, RegulationType>();
	}

	public void addCatalyst(Entity e) {
		catalyst.add(e);
	}

	public void addInput(Entity e) {
		inputs.add(e);
	}

	public void addOutput(Entity e) {
		outputs.add(e);
	}

	public void addRegulator(Entity e, RegulationType type) {
		regulator.put(e, type);
	}

	public void addRegulator(Entity e, String type) {
		regulator.put(e, RegulationType.valueOf(type));
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

		// Calculate inputs, catalysts & regulators
		for (Entity ein : inputs)
			ein.calc(doneEntities);

		for (Entity ecat : catalyst)
			ecat.calc(doneEntities);

		for (Entity ereg : regulator.keySet())
			ereg.calc(doneEntities);

		// Calculate aggregated input
		double in = Double.POSITIVE_INFINITY;
		for (Entity ein : inputs)
			if (ein.hasOutput()) in = Math.min(in, ein.getOutput());

		// Apply 'catalyst'
		double cat = 1.0; // Neutral by default
		for (Entity ecat : catalyst) {
			if (ecat.hasOutput()) {
				double dcat = ecat.getOutput();
				double sigm = 2.0 / (1.0 + Math.exp(-dcat));
				cat *= sigm;
			}
		}

		// Apply 'regulation'
		double reg = 1.0; // Neutral by default
		for (Entity ereg : regulator.keySet()) {
			if (ereg.hasOutput()) {
				double dcat = ereg.getOutput();
				double sigm = 1.0;

				RegulationType regType = regulator.get(ereg);

				switch (regType) {

				case PositiveRegulation:
					sigm = 2.0 / (1.0 + Math.exp(-dcat));
					break;

				case NegativeRegulation:
					sigm = 2.0 / (1.0 + Math.exp(dcat));
					break;

				case Requirement:
					sigm = 1.0 / (1.0 + Math.exp(dcat));
					break;
				}

				reg *= sigm; // Summarize weight
			}
		}

		// Nothing in input? => Cannot calculate output
		if (Double.isInfinite(in)) output = Double.NaN;
		else {
			output = in * cat * reg;
			Gpr.debug("REACTION:\t" + this.getName() + "\t" + in + "\t" + cat + "\t" + reg + "\t=>\t" + output);
		}

		return output;
	}

	public HashSet<Entity> getInputs() {
		return inputs;
	}

	public HashSet<Entity> getOutputs() {
		return outputs;
	}

	public HashMap<Entity, RegulationType> getRegulator() {
		return regulator;
	}

	@Override
	public String toString() {
		return toString(0, new HashSet<Entity>());
	}

	@Override
	public String toString(int tabs, HashSet<Entity> done) {
		StringBuilder sb = new StringBuilder();
		sb.append(Gpr.tabs(tabs) + getClass().getSimpleName() + "[" + id + "]: " + name + "\n");

		if (!inputs.isEmpty()) {
			sb.append(Gpr.tabs(tabs + 1) + "Inputs:\n");
			for (Entity e : inputs) {
				if (done.contains(e)) sb.append(Gpr.tabs(tabs + 2) + e.toStringSimple() + "\n");
				else {
					done.add(e);
					sb.append(e.toString(tabs + 2, done) + "\n");
				}
			}
		}

		if (!outputs.isEmpty()) {
			sb.append(Gpr.tabs(tabs + 1) + "Outputs:\n");
			for (Entity e : outputs) {
				if (done.contains(e)) sb.append(Gpr.tabs(tabs + 2) + e.toStringSimple() + "\n");
				else {
					done.add(e);
					sb.append(e.toString(tabs + 2, done) + "\n");
				}
			}
		}

		if (!catalyst.isEmpty()) {
			sb.append(Gpr.tabs(tabs + 1) + "Catalysts:\n");
			for (Entity e : catalyst) {
				if (done.contains(e)) sb.append(Gpr.tabs(tabs + 2) + e.toStringSimple() + "\n");
				else {
					done.add(e);
					sb.append(e.toString(tabs + 2, done) + "\n");
				}
			}
		}

		if (!regulator.isEmpty()) {
			sb.append(Gpr.tabs(tabs + 1) + "Regulator:\n");
			for (Entity e : regulator.keySet()) {
				if (done.contains(e)) sb.append(Gpr.tabs(tabs + 2) + e.toStringSimple() + "(" + regulator.get(e) + ")\t" + "\n");
				else {
					done.add(e);
					sb.append(e.toString(tabs + 2, done) + "(" + regulator.get(e) + ")\t" + "\n");
				}
			}
		}

		return sb.toString();
	}

}
