package ca.mcgill.mcb.pcingola.snpEffect.testCases;

import java.util.Random;

import junit.framework.Assert;
import junit.framework.TestCase;
import ca.mcgill.mcb.pcingola.gsa.GenePvalueList;
import ca.mcgill.mcb.pcingola.gsa.GenePvalueList.PvalueSummary;
import flanagan.analysis.Stat;

/**
 * GenePvalueList statistics test case
 * 
 * @author pcingola
 */
public class TestGenePvalueList extends TestCase {

	/**
	 * Combined p-value : MIN
	 */
	public void test_01() {
		double pvals[] = { 0.01, 0.2, 0.3 };

		// Create p values
		GenePvalueList gpl = new GenePvalueList();
		for (double pval : pvals)
			gpl.add(pval);

		// Check pvalues
		double pvalue = gpl.pValue(PvalueSummary.MIN);
		Assert.assertEquals(0.01, pvalue);
	}

	/**
	 * Combined p-value : AVG
	 */
	public void test_02() {
		double pvals[] = { 0.01, 0.2, 0.3 };

		// Create p values
		GenePvalueList gpl = new GenePvalueList();
		for (double pval : pvals)
			gpl.add(pval);

		// Check pvalues
		double pvalue = gpl.pValue(PvalueSummary.AVG);
		Assert.assertEquals(0.17, pvalue);
	}

	/**
	 * Combined p-value : AVG10
	 */
	public void test_03() {
		double pvals[] = { 0.01, 0.9, 0.2, 0.9, 0.3, 0.9, 0.01, 0.9, 0.2, 0.9, 0.3, 0.9, 0.01, 0.9, 0.2, 0.9, 0.3, 0.9, 0.17, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9 };

		// Create p values
		GenePvalueList gpl = new GenePvalueList();
		for (double pval : pvals)
			gpl.add(pval);

		// Check pvalues
		double pvalue = gpl.pValue(PvalueSummary.AVG10);
		Assert.assertEquals(0.17, pvalue);
	}

	/**
	 * Complementary CDF for Chi^2 distribution
	 * 
	 */
	public void test_04() {
		Random rand = new Random(20130609);

		// Create many random tests
		for (int i = 0; i < 100000; i++) {
			// Random degrees of freedom
			int degOfFreedom = rand.nextInt(20) + 1;

			// Random chi^2
			double chi2 = 0;
			for (int j = 0; j < degOfFreedom; j++) {
				double z = rand.nextGaussian();
				chi2 += z * z;
			}

			// Calculate complementary probabilities
			double pval = GenePvalueList.chiSquareCDFComplementary(chi2, degOfFreedom);
			double prob = Stat.chiSquareCDF(chi2, degOfFreedom);

			// Assert that statistics add to 1.0
			Assert.assertEquals(1.0, pval + prob);
		}
	}

	/**
	 * Combined p-value : FISHER_CHI_SQUARE
	 */
	public void test_05() {
		double pvals[] = { 0.01, 0.2, 0.3 };

		// Create p values
		GenePvalueList gpl = new GenePvalueList();
		for (double pval : pvals)
			gpl.add(pval);

		// Check pvalues
		double pvalue = gpl.pValue(PvalueSummary.FISHER_CHI_SQUARE);
		Assert.assertEquals(0.02156175132483462, pvalue);
	}

	/**
	 * Combined p-value : Z_SCORES
	 */
	public void test_06() {
		double pvals[] = { 0.01, 0.2, 0.3 };

		// Create p values
		GenePvalueList gpl = new GenePvalueList();
		for (double pval : pvals)
			gpl.add(pval);

		// Check pvalues
		double pvalue = gpl.pValue(PvalueSummary.Z_SCORES);
		Assert.assertEquals(0.01651203252368999, pvalue);
	}

}
