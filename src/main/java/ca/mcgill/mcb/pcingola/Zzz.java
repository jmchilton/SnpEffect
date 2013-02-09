package ca.mcgill.mcb.pcingola;

import ca.mcgill.mcb.pcingola.codons.CodonTable;
import ca.mcgill.mcb.pcingola.codons.CodonTables;

/**
 * Simple test program
 * @author pcingola
 */
public class Zzz {

	public static void main(String[] args) {
		CodonTable codonTable = CodonTables.getInstance().getTable(CodonTables.STANDARD_TABLE_NAME);

		System.out.println(codonTable);
		System.out.println(codonTable.aaThreeLetterCode("*"));
	}
}
