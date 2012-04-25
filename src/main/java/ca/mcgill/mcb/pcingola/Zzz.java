package ca.mcgill.mcb.pcingola;

import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;

public class Zzz {

	public static void main(String[] args) {
		String vcfFileName = Gpr.HOME + "/snpEFf/empty.vcf";
		VcfFileIterator vcfFile = new VcfFileIterator(vcfFileName);
		for (VcfEntry ve : vcfFile) {
			System.out.println(ve.isVariant() + "\t" + ve);
		}

	}
}
