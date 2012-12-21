package ca.mcgill.mcb.pcingola.snpEffect.testCases;

import java.util.Random;

import junit.framework.Assert;
import junit.framework.TestCase;
import ca.mcgill.mcb.pcingola.genotypes.GenotypeVector;

/**
 * Test cases for GenotypeVector class
 * 
 * @author pcingola
 */
public class TestCasesGenotypeVector extends TestCase {

	public void test_01() {
		// Show masks (just to check they are OK)
		for (byte m : GenotypeVector.mask)
			System.out.println("Mask          :" + m + "\t" + Integer.toBinaryString(m & 0xff));

		for (byte m : GenotypeVector.reverseMask)
			System.out.println("Reverse Mask  :" + m + "\t" + Integer.toBinaryString(m & 0xff));

		for (int code = 0; code < 4; code++) {
			GenotypeVector gv = new GenotypeVector(2);

			for (int i = 0; i < 4; i++)
				gv.set(i, code);

			for (int i = 0; i < 4; i++)
				Assert.assertEquals(code, gv.get(i));
		}
	}

	public void test_02() {
		Random rand = new Random(20121221);
		GenotypeVector gv = new GenotypeVector(1000);

		// Create random codes
		int codes[] = new int[gv.size()];
		for (int i = 0; i < gv.size(); i++) {
			int code = rand.nextInt(4);
			codes[i] = code;
			gv.set(i, code);
		}

		// Check that codes are stored OK
		for (int i = 0; i < gv.size(); i++) {
			Assert.assertEquals(codes[i], gv.get(i));
		}
	}
}
