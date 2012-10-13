package ca.mcgill.mcb.pcingola;

import ca.mcgill.mcb.pcingola.fileIterator.FastaFileIterator;
import ca.mcgill.mcb.pcingola.util.Gpr;

public class Zzz {

	public static void main(String[] args) {
		String file = Gpr.HOME + "/snpEff/data/genomes/GRCm38.68.fa";

		System.out.println("\t\tReading file '" + file);
		FastaFileIterator ffi = new FastaFileIterator(file);
		for (String seq : ffi) {
			String chromo = ffi.getName();
			System.out.println("\t\tReading sequence '" + chromo + "', length: " + seq.length());
		}

	}
}
