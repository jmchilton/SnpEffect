package ca.mcgill.mcb.pcingola.snpEffect.testCases;

import java.util.Random;

import junit.framework.Assert;
import junit.framework.TestCase;
import ca.mcgill.mcb.pcingola.gsa.PvalueList;
import ca.mcgill.mcb.pcingola.gsa.PvalueList.PvalueSummary;
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
		PvalueList gpl = new PvalueList();
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
		PvalueList gpl = new PvalueList();
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
		PvalueList gpl = new PvalueList();
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
			double pval = PvalueList.chiSquareCDFComplementary(chi2, degOfFreedom);
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
		PvalueList gpl = new PvalueList();
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
		PvalueList gpl = new PvalueList();
		for (double pval : pvals)
			gpl.add(pval);

		// Check pvalues
		double pvalue = gpl.pValue(PvalueSummary.Z_SCORES);
		Assert.assertEquals(0.01651203252368999, pvalue);
	}

	/**
	 * Combined p-value : FDR
	 */
	public void test_07() {
		double pvals[] = { 2.354054e-07, 2.101590e-05, 2.576842e-05, 9.814783e-05, 1.052610e-04, 1.241481e-04, 1.325988e-04, 1.568503e-04, 2.254557e-04, 3.795380e-04, 6.114943e-04, 1.613954e-03, 3.302430e-03, 3.538342e-03, 5.236997e-03, 6.831909e-03, 7.059226e-03, 8.805129e-03, 9.401040e-03, 1.129798e-02, 2.115017e-02, 4.922736e-02, 6.053298e-02, 6.262239e-02, 7.395153e-02, 8.281103e-02, 8.633331e-02, 1.190654e-01, 1.890796e-01, 2.058494e-01, 2.209214e-01, 2.856000e-01, 3.048895e-01, 4.660682e-01, 4.830809e-01, 4.921755e-01, 5.319453e-01, 5.751550e-01, 5.783195e-01, 6.185894e-01, 6.363620e-01, 6.448587e-01, 6.558414e-01, 6.885884e-01, 7.189864e-01, 8.179539e-01, 8.274487e-01, 8.971300e-01, 9.118680e-01, 9.437890e-01 };

		// Create p values
		PvalueList gpl = new PvalueList();
		for (double pval : pvals)
			gpl.add(pval);

		// Check pvalues
		double pvalue = gpl.pValue(PvalueSummary.FDR);
		Assert.assertEquals(0.028244949999999998, pvalue);
	}

}
