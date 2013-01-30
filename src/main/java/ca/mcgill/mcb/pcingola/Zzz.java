package ca.mcgill.mcb.pcingola;

import ca.mcgill.mcb.pcingola.fileIterator.FastaFileIterator;

/**
 * Simple test program
 * @author pcingola
 */
public class Zzz {

	public static void main(String[] args) {
		String fastaFileName = "/Users/pablocingolani/snpEff/z.fa";
		FastaFileIterator ffi = new FastaFileIterator(fastaFileName);
		for (String seq : ffi) {
			System.out.println("SeqName: " + ffi.getName() + "\tSize: " + seq.length());
		}

	}
}
