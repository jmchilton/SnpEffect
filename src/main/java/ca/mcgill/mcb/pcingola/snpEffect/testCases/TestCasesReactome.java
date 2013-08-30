package ca.mcgill.mcb.pcingola.snpEffect.testCases;

import java.util.HashSet;

import junit.framework.Assert;
import junit.framework.TestCase;
import ca.mcgill.mcb.pcingola.reactome.Entity;
import ca.mcgill.mcb.pcingola.reactome.Entity.TransferFunction;
import ca.mcgill.mcb.pcingola.reactome.events.Reaction;
import ca.mcgill.mcb.pcingola.reactome.events.Reaction.RegulationType;

/**
 * Test Reactome circuits
 * 
 * @author pcingola
 */
public class TestCasesReactome extends TestCase {

	public static boolean debug = false;
	public static boolean verbose = true;
	public static int SHOW_EVERY = 10;

	public TestCasesReactome() {
		super();
	}

	/**
	 * Reaction with two molecules
	 */
	public void test_01() {
		int id = 1;
		Entity e1 = new Entity(id++, "input_1");
		Entity e2 = new Entity(id++, "input_2");

		Reaction r = new Reaction(id++, "reaction_1");
		r.addInput(e1);
		r.addInput(e2);

		e1.setFixedOutput(3.2);
		e2.setFixedOutput(2.1);

		Entity.debug = true;

		Entity.TRANSFER_FUNCTION = TransferFunction.LINEAR;
		double out = r.calc(new HashSet<Entity>());
		System.out.println("Out: " + out);
		Assert.assertEquals(2.1, out);

		Entity.TRANSFER_FUNCTION = TransferFunction.TANH;
		out = r.calc(new HashSet<Entity>());
		System.out.println("Out: " + out);
		Assert.assertEquals(0.9999501692221211, out);

	}

	/**
	 * Reaction with a catalyst
	 */
	public void test_02() {
		int id = 1;
		Entity e1 = new Entity(id++, "input_1");
		Entity e2 = new Entity(id++, "input_2");
		Entity cat = new Entity(id++, "catalyst");

		Reaction r = new Reaction(id++, "reaction_1");
		r.addInput(e1);
		r.addInput(e2);
		r.addCatalyst(cat);

		e1.setFixedOutput(3.2);
		e2.setFixedOutput(2.1);
		cat.setWeight(-1.0);

		Entity.TRANSFER_FUNCTION = TransferFunction.LINEAR;

		double out = r.calc(new HashSet<Entity>());
		System.out.println("Out: " + out);
		Assert.assertEquals(1.1295539697539796, out);

		cat.setWeight(2.0);
		out = r.calc(new HashSet<Entity>());
		System.out.println("Out: " + out);
		Assert.assertEquals(3.699347727507106, out);

		Entity.TRANSFER_FUNCTION = TransferFunction.TANH;
		e1.setFixedOutput(0.9);
		e2.setFixedOutput(0.5);
		cat.setWeight(-1.0);
		out = r.calc(new HashSet<Entity>());
		System.out.println("Out: " + out);
		Assert.assertEquals(0.6369566834708253, out);
	}

	/**
	 * Reaction with positive regulation
	 */
	public void test_03() {
		int id = 1;
		Entity e1 = new Entity(id++, "input_1");
		Entity e2 = new Entity(id++, "input_2");
		Entity reg = new Entity(id++, "catalyst");

		Reaction r = new Reaction(id++, "reaction_1");
		r.addInput(e1);
		r.addInput(e2);
		r.addRegulator(reg, RegulationType.PositiveRegulation);

		e1.setFixedOutput(3.2);
		e2.setFixedOutput(2.1);
		reg.setFixedOutput(-1.0);

		Entity.TRANSFER_FUNCTION = TransferFunction.LINEAR;

		double out = r.calc(new HashSet<Entity>());
		System.out.println("Out: " + out);
		Assert.assertEquals(2.66477698487699, out);

		reg.setFixedOutput(2.0);
		out = r.calc(new HashSet<Entity>());
		System.out.println("Out: " + out);
		Assert.assertEquals(3.949673863753553, out);
	}

	/**
	 * Reaction with negative regulation
	 */
	public void test_04() {
		int id = 1;
		Entity e1 = new Entity(id++, "input_1");
		Entity e2 = new Entity(id++, "input_2");
		Entity reg = new Entity(id++, "catalyst");

		Reaction r = new Reaction(id++, "reaction_1");
		r.addInput(e1);
		r.addInput(e2);
		r.addRegulator(reg, RegulationType.NegativeRegulation);

		e1.setFixedOutput(3.2);
		e2.setFixedOutput(2.1);
		reg.setFixedOutput(-1.0);

		Entity.TRANSFER_FUNCTION = TransferFunction.LINEAR;

		double out = r.calc(new HashSet<Entity>());
		System.out.println("Out: " + out);
		Assert.assertEquals(1.5352230151230104, out);

		reg.setFixedOutput(2.0);
		out = r.calc(new HashSet<Entity>());
		System.out.println("Out: " + out);
		Assert.assertEquals(0.2503261362464472, out);
	}

	/**
	 * Reaction with requirement
	 */
	public void test_05() {
		int id = 1;
		Entity e1 = new Entity(id++, "input_1");
		Entity e2 = new Entity(id++, "input_2");
		Entity reg = new Entity(id++, "catalyst");

		Reaction r = new Reaction(id++, "reaction_1");
		r.addInput(e1);
		r.addInput(e2);
		r.addRegulator(reg, RegulationType.Requirement);

		e1.setFixedOutput(3.2);
		e2.setFixedOutput(2.1);
		reg.setFixedOutput(-1.0);

		Entity.TRANSFER_FUNCTION = TransferFunction.LINEAR;

		double out = r.calc(new HashSet<Entity>());
		System.out.println("Out: " + out);
		Assert.assertEquals(0.5647769848769898, out);

		reg.setFixedOutput(2.0);
		out = r.calc(new HashSet<Entity>());
		System.out.println("Out: " + out);
		Assert.assertEquals(1.849673863753553, out);
	}

}
