package ca.mcgill.mcb.pcingola.reactome.events;

import java.util.Collection;
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

	public static final double MAX_WEIGHT_SUM = 1.1;

	protected HashMap<Entity, Double> inputs; // Input -> Weight
	protected HashSet<Entity> outputs;
	protected HashSet<Entity> catalyst;
	protected HashMap<Entity, RegulationType> regulator;

	public Reaction(int id, String name) {
		super(id, name);
		inputs = new HashMap<Entity, Double>();
		outputs = new HashSet<Entity>();
		catalyst = new HashSet<Entity>();
		regulator = new HashMap<Entity, RegulationType>();
	}

	public void addCatalyst(Entity e) {
		catalyst.add(e);
	}

	public void addInput(Entity e) {
		if (e == null) return;

		if (!inputs.containsKey(e)) inputs.put(e, 1.0);
		else inputs.put(e, inputs.get(e) + 1.0);
	}

	public void addOutput(Entity e) {
		outputs.add(e);

		// Add input
		if (e.isReaction()) ((Reaction) e).addInput(this);
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

		if (!Double.isNaN(fixedOutput)) output = fixedOutput;
		else {
			if (doneEntities.contains(this)) return output; // Make sure we don't calculate twice
			doneEntities.add(this); // Keep 'entities' set up to date

			// Calculate inputs
			for (Entity ein : getInputs())
				ein.calc(doneEntities);

			// Calculate catalysts 
			for (Entity ecat : catalyst)
				ecat.calc(doneEntities);

			// Calculate regulators
			for (Entity ereg : regulator.keySet())
				ereg.calc(doneEntities);

			// Calculate aggregated input
			double in = 0;
			switch (TRANSFER_FUNCTION) {
			case LINEAR:
				// Min of all inputs
				in = Double.POSITIVE_INFINITY;
				for (Entity ein : getInputs())
					if (ein.hasOutput()) in = Math.min(in, ein.getOutput());
				break;

			default:
				// Weighted input
				in = 0;
				for (Entity ein : getInputs())
					if (ein.hasOutput()) in += ein.getOutput() * inputs.get(ein);
				break;
			}

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
						sigm = 1 + 1.0 / (1.0 + Math.exp(-dcat));
						break;

					case NegativeRegulation:
						sigm = 1 - 1.0 / (1.0 + Math.exp(-dcat));
						break;

					case Requirement:
						sigm = 1.0 / (1.0 + Math.exp(-dcat));
						break;
					}

					reg *= sigm; // Summarize weight
				}
			}

			// Nothing in input? => Cannot calculate output
			if (Double.isInfinite(in)) output = Double.NaN;
			else output = transferFunction(in, cat, reg);
		}

		if (debug) System.out.println(output + "\tfixed:" + isFixed() + "\tid:" + id + "\ttype:" + getClass().getSimpleName() + "\tname:" + name);
		return output;
	}

	/**
	 * Transfer function
	 * @param x
	 * @return
	 */
	protected double transferFunction(double x, double cat, double reg) {
		switch (TRANSFER_FUNCTION) {
		case SIGM_PLUS_MINUS:
			return 2.0 * (cat * reg) / (1.0 + Math.exp(-x)) - 1.0;
		case SIGM:
			return (cat * reg) / (1.0 + Math.exp(-x));
		case LINEAR:
		default:
			throw new RuntimeException("Unimplemented transfer function: " + TRANSFER_FUNCTION);
		}
	}

	public HashSet<Entity> getCatalyst() {
		return catalyst;
	}

	public Collection<Entity> getInputs() {
		return inputs.keySet();
	}

	public HashSet<Entity> getOutputs() {
		return outputs;
	}

	public HashMap<Entity, RegulationType> getRegulator() {
		return regulator;
	}

	@Override
	public boolean isReaction() {
		return true;
	}

	public void scaleWeights() {
		// One input? Nothing to do
		if (inputs.size() <= 1) return;

		// Sum
		double sum = 0;
		for (Double in : inputs.values())
			sum += Math.abs(in);

		// Is it scaled? Nothing to do
		if (sum < MAX_WEIGHT_SUM) return;

		// New hash with scaled weights
		HashMap<Entity, Double> newInputs = new HashMap<Entity, Double>();
		for (Entity e : inputs.keySet()) {
			double weight = inputs.get(e) / sum;
			newInputs.put(e, weight);
		}

		// Replace hash
		inputs = newInputs;

	}

	@Override
	public String toString() {
		return toString(0, new HashSet<Entity>());
	}

	@Override
	public String toString(int tabs, HashSet<Entity> done) {
		done.add(this);

		StringBuilder sb = new StringBuilder();
		sb.append(Gpr.tabs(tabs) + getClass().getSimpleName() + "[" + id + "]: " + name + "\n");

		if (!inputs.isEmpty()) {
			sb.append(Gpr.tabs(tabs + 1) + "Inputs:\n");
			for (Entity e : getInputs()) {
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
