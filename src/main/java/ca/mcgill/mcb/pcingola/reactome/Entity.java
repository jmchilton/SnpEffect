package ca.mcgill.mcb.pcingola.reactome;

import java.util.Collection;
import java.util.HashSet;

import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * A reactome basic entity (e.g. anything in reactome database derives fro this object)
 * 
 * @author pcingola
 *
 */
public class Entity {

	public enum TransferFunction {
		LINEAR, TANH, SIGM
	};

	public static boolean debug = false;
	public static TransferFunction TRANSFER_FUNCTION = TransferFunction.TANH;

	protected int id; // Entity ID
	protected String name; // Entity Name
	protected Compartment compartment; // In which cell's compartment is this entity located?
	protected double output; // Entity output value
	protected double weight; // Weight applied to this entity (NaN means not available)
	protected double fixedOutput; // Fixed output value (external information)
	protected int countWeights;
	protected HashSet<String> geneIds; // All gene IDs related to this entity

	public static final double MAX_OUTPUT = 10;
	public static final double MIN_OUTPUT = -10;

	public Entity(int id, String name) {
		this.id = id;
		this.name = name;
		reset();
	}

	/**
	 * Add a geneId
	 * @param geneId
	 */
	public void addGeneId(String geneId) {
		if (geneIds == null) geneIds = new HashSet<String>();
		geneIds.add(geneId);
	}

	public void addWeight(double weight) {
		if (Double.isNaN(weight)) return; // Nothing to do
		else if (countWeights == 0) this.weight = weight; // In case previous value was NAN
		else this.weight += weight;
		countWeights++;
	}

	/**
	 * Calculate entities.
	 * Make sure we don't calculate twice (keep 'doneEntities' set up to date)
	 * 
	 * @param doneEntities
	 * @return
	 */
	public double calc(HashSet<Entity> doneEntities) {

		if (!Double.isNaN(fixedOutput)) {
			output = fixedOutput;
		} else {
			if (doneEntities.contains(this)) return output; // Make sure we don't calculate twice
			doneEntities.add(this); // Keep 'entities' set up to date

			output = getWeight(); // Calculate output
		}

		if (debug) System.out.println(output + "\tfixed:" + isFixed() + "\tid:" + id + "\ttype:" + getClass().getSimpleName() + "\tname:" + name);
		return output;
	}

	public Compartment getCompartment() {
		return compartment;
	}

	public Collection<String> getGeneIds() {
		return geneIds;
	}

	public int getId() {
		return id;
	}

	public String getName() {
		return name;
	}

	public double getOutput() {
		return output;
	}

	public double getWeight() {
		if (countWeights > 1) return weight / Math.sqrt(countWeights);
		return weight;
	}

	public boolean hasOutput() {
		return !Double.isNaN(output);
	}

	public boolean isFixed() {
		return !Double.isNaN(fixedOutput);
	}

	public boolean isReaction() {
		return false;
	}

	public void reset() {
		output = Double.NaN; // Entity output value
		weight = Double.NaN; // Weight applied to this entity (NaN means not available)
		fixedOutput = Double.NaN; // Fixed output value (external information)
		countWeights = 0;

		output = weight = 0;
	}

	public void resetWeight() {
		weight = Double.NaN;
		countWeights = 0;
	}

	public void setCompartment(Compartment compartment) {
		this.compartment = compartment;
	}

	public void setFixedOutput(double fixedOutput) {
		this.fixedOutput = fixedOutput;
	}

	public void setWeight(double weight) {
		this.weight = weight;
		countWeights = 1;
	}

	@Override
	public String toString() {
		return toString(0, new HashSet<Entity>());
	}

	public String toString(int tabs, HashSet<Entity> done) {
		done.add(this);
		return Gpr.tabs(tabs) + getClass().getSimpleName() + "[" + id + "]: " + name;
	}

	public String toStringSimple() {
		return getClass().getSimpleName() + "[" + id + "]: " + name;
	}

	/**
	 * Transfer function
	 * @param x
	 * @return
	 */
	protected double transferFunction(double x) {
		switch (TRANSFER_FUNCTION) {
		case LINEAR:
			return x;
		case SIGM:
			return 1.0 / (1.0 + Math.exp(-x));
		case TANH:
			return Math.tanh(x);
		default:
			throw new RuntimeException("Unimplemented transfer function: " + TRANSFER_FUNCTION);
		}
	}

}
