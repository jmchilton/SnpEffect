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

	protected int id; // Entity ID
	protected String name; // Entity Name
	protected Compartment compartment; // In which cell's compartment is this entity located?
	protected double output = Double.NaN; // Entity output value
	protected double weight = Double.NaN; // Weight applied to this entity (NaN means not available)
	protected int countWeights = 0;
	protected HashSet<String> geneIds; // All gene IDs related to this entity

	public Entity(int id, String name) {
		this.id = id;
		this.name = name;
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
		if (doneEntities.contains(this)) return output; // Make sure we don't calculate twice

		doneEntities.add(this); // Keep 'entities' set up to date
		output = getWeight(); // Calculate output

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
		if (countWeights > 1) {
			weight = weight / countWeights;
			countWeights = 1;
		}
		return weight;
	}

	public boolean hasOutput() {
		return !Double.isNaN(output);
	}

	public boolean isEvent() {
		return false;
	}

	public void setCompartment(Compartment compartment) {
		this.compartment = compartment;
	}

	public void setWeight(double weight) {
		this.weight = weight;
		countWeights = 1;
	}

	@Override
	public String toString() {
		return getClass().getSimpleName() + "[" + id + "]: " + name;
	}

	public String toString(int tabs, HashSet<Entity> done) {
		return Gpr.tabs(tabs) + getClass().getSimpleName() + "[" + id + "]: " + name;
	}

	public String toStringSimple() {
		return getClass().getSimpleName() + "[" + id + "]: " + name;
	}

}
