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

	public void test_04() {
		Random rand = new Random(20130609);

		for (int i = 0; i < 100000; i++) {
			int k = rand.nextInt(20) + 1;

			double chi2 = 0;
			for (int j = 0; j < k; j++) {
				double z = rand.nextGaussian();
				chi2 += z * z;
			}

			// Calculate complementary probabilities
			double pval = GenePvalueList.chiSquareCDFComplementary(chi2, k);
			double prob = Stat.chiSquareCDF(chi2, k);

			// Gpr.debug("Chi2 : " + chi2 + "\tdf: " + k + "\t" + pval + "\t" + prob + "\t" + (pval + prob));
			Assert.assertEquals(1.0, pval + prob);
		}
	}

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

}
