package ca.mcgill.mcb.pcingola.snpEffect.testCases;

import java.io.IOException;
import java.util.Random;

import junit.framework.Assert;
import junit.framework.TestCase;
import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.vcf.FileIndexChrPos;
import ca.mcgill.mcb.pcingola.vcf.FileIndexChrPos.LineAndPos;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;

/**
 * Test cases for file index (chr:pos index on files)
 * 
 * @author pcingola
 */
public class TestCasesFileIndexChrPos extends TestCase {

	boolean verbose = true;
	boolean debug = false;

	/**
	 * Test : Find beginning of a chromosome
	 * @throws IOException
	 */
	public void test_01() throws IOException {
		String vcf = "tests/test.chr1.vcf";

		System.out.println("Indexing file '" + vcf + "'");
		FileIndexChrPos idx = new FileIndexChrPos(vcf);
		idx.setVerbose(verbose);
		idx.open();
		idx.index();

		long pos = idx.getStart("1");
		System.out.println("\tChr 1 start: " + pos);
		Assert.assertEquals(82703, pos);

		idx.close();
	}

	/**
	 * Test : Find a line
	 * @throws IOException
	 */
	public void test_02() throws IOException {
		String vcf = "tests/test.chr1.vcf";
		String line = "1	861275	.	C	T	764.18	PASS	AC=1;AF=0.00061;AN=1644;DS;set=Intersection";

		System.out.println("Indexing file '" + vcf + "'");
		FileIndexChrPos idx = new FileIndexChrPos(vcf);
		idx.setVerbose(verbose);
		idx.open();
		idx.index();

		long pos = idx.getStart("1");
		LineAndPos lp = idx.getLine(pos);
		System.out.println("\tChr 1 start: " + pos + "\tLine: '" + lp.line + "'");
		Assert.assertEquals(line, lp.line);
		idx.close();
	}

	/**
	 * Test : Find a chr:pos
	 * @throws IOException
	 */
	public void test_03() throws IOException {
		String vcf = "tests/test.chr1.vcf";

		System.out.println("Indexing file '" + vcf + "'");
		FileIndexChrPos idx = new FileIndexChrPos(vcf);
		idx.setVerbose(verbose);
		idx.open();
		idx.index();

		int chrPos = 861275 - 1; // Zero based coordinate of position in first VCF line

		// Beginning of line
		long pos = idx.find("1", chrPos, false);
		Gpr.debug("Pos: " + pos);
		Assert.assertEquals(82703, pos);

		// End of line
		pos = idx.find("1", chrPos, true);
		Gpr.debug("Pos: " + pos);
		Assert.assertEquals(82774, pos);

		idx.close();
	}

	/**
	 * Test : Find a chr:pos
	 * @throws IOException
	 */
	public void test_04() throws IOException {
		String vcf = "tests/test.chr1.vcf";

		System.out.println("Indexing file '" + vcf + "'");
		FileIndexChrPos idx = new FileIndexChrPos(vcf);
		idx.setVerbose(verbose);
		idx.setDebug(debug);
		idx.open();
		idx.index();

		// We'll try to find this chr:pos = 1:1019717 
		int chrPos = 1019717 - 1; // Zero based coordinate of VCF line

		// Find chr:pos
		long pos = idx.find("1", chrPos, false);
		LineAndPos lp = idx.getLine(pos);
		int chrPosLp = idx.pos(lp.line);
		Assert.assertEquals(chrPos, chrPosLp);
		Assert.assertEquals(129869, pos);

		idx.close();
	}

	/**
	 * Test : Find a chr:pos that does not exists
	 * @throws IOException
	 */
	public void test_05() throws IOException {
		String vcf = "tests/test.chr1.vcf";

		System.out.println("Indexing file '" + vcf + "'");
		FileIndexChrPos idx = new FileIndexChrPos(vcf);
		idx.setVerbose(verbose);
		idx.setDebug(debug);
		idx.open();
		idx.index();

		// We'll try to find this chr:pos = 1:1019716 (the coordinate that is in the VCF file is 1:1019717) 
		int chrPosReal = 1019717 - 1; // Zero based coordinate of VCF line
		int chrPos = chrPosReal - 1; // Zero based coordinate of VCF line

		// Find chr:pos
		long pos = idx.find("1", chrPos, false);
		LineAndPos lp = idx.getLine(pos);
		int chrPosLp = idx.pos(lp.line);
		Gpr.debug(lp.line);
		Assert.assertEquals(chrPosReal, chrPosLp);

		idx.close();
	}

	/**
	 * Test : Find a chr:pos that does not exists
	 * @throws IOException
	 */
	public void test_06() throws IOException {
		String vcf = "tests/test.chr1.vcf";

		System.out.println("Indexing file '" + vcf + "'");
		FileIndexChrPos idx = new FileIndexChrPos(vcf);
		idx.setVerbose(verbose);
		idx.setDebug(debug);
		idx.open();
		idx.index();

		// We'll try to find this chr:pos = 1:1019716 (the coordinate that is in the VCF file is 1:1019717) 
		int chrPosReal = 1019718 - 1; // Zero based coordinate of VCF line
		int chrPos = chrPosReal + 1; // Zero based coordinate of VCF line

		// Find chr:pos
		long pos = idx.find("1", chrPos, false);
		LineAndPos lp = idx.getLine(pos);
		int chrPosLp = idx.pos(lp.line);
		Assert.assertEquals(1019733 - 1, chrPosLp); // We expect to find the next coordinate in VCF file (zero-based)

		idx.close();
	}

	/**
	 * Test : Find a chr:pos that does not exists
	 * @throws IOException
	 */
	public void test_07() throws IOException {
		String vcf = "tests/test.chr1.vcf";

		System.out.println("Indexing file '" + vcf + "'");
		FileIndexChrPos idx = new FileIndexChrPos(vcf);
		idx.setVerbose(verbose);
		idx.setDebug(debug);
		idx.open();
		idx.index();

		// We'll try to find this chr:pos = 1:1019716 (the coordinate that is in the VCF file is 1:1019717) 
		int chrPosReal = 865488 - 1; // Zero based coordinate of VCF line
		int chrPos = 861316 - 1; // Zero based coordinate of VCF line

		// Find chr:pos
		long pos = idx.find("1", chrPos, false);
		LineAndPos lp = idx.getLine(pos);
		int chrPosLp = idx.pos(lp.line);
		Assert.assertEquals(chrPosReal, chrPosLp); // We expect to find the next coordinate in VCF file (zero-based)

		idx.close();
	}

	/**
	 * Test : Find a chr:pos that does not exists
	 * @throws IOException
	 */
	public void test_10() throws IOException {
		String vcfFileName = "tests/test.chr1.vcf";

		Random random = new Random(20130217);

		// Index file
		System.out.println("Indexing file '" + vcfFileName + "'");
		FileIndexChrPos idx = new FileIndexChrPos(vcfFileName);
		idx.setVerbose(verbose);
		idx.setDebug(debug);
		idx.open();
		idx.index();

		// Iterate over vcf file
		VcfFileIterator vcf = new VcfFileIterator(vcfFileName);
		int chrPosPrev = 0;
		for (VcfEntry ve : vcf) {
			int chrPos = ve.getStart();
			if (chrPosPrev == 0) chrPosPrev = chrPos - 100;

			// Only perform some tests otherwise it's too long
			if (random.nextInt(1000) < 2) {
				Gpr.debug("Find: " + chrPosPrev + " - " + chrPos);

				// Find all positions from previous to current
				int step = Math.max((chrPos - chrPosPrev) / 10, 1);
				for (int cp = chrPosPrev; cp <= chrPos; cp += step) {
					// Find chr:pos
					long pos = idx.find("1", cp, false);
					LineAndPos lp = idx.getLine(pos);
					int chrPosLp = idx.pos(lp.line);
					if (debug) Gpr.debug("Find: " + cp + "\t" + lp.line);
					Assert.assertEquals(chrPos, chrPosLp); // We expect to find the next coordinate in VCF file (zero-based)
				}
			}

			chrPosPrev = chrPos + 1;
		}

		idx.close();
	}
}
