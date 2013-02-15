package ca.mcgill.mcb.pcingola;

import java.util.List;

import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.vcf.PedigreeEnrty;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;

/**
 * Simple test program
 * @author pcingola
 */
public class Zzz {

	public static void main(String[] args) {
		String file = Gpr.HOME + "/snpEff/test.cancer.snp.vcf";
		VcfFileIterator vcf = new VcfFileIterator(file);

		for (VcfEntry ve : vcf) {
			if (vcf.isHeadeSection()) {
				Gpr.debug("HEADER");
				List<PedigreeEnrty> pedigre = vcf.getVcfHeader().getPedigree();
				for (PedigreeEnrty pe : pedigre)
					Gpr.debug("\t" + pe);
			}

			Gpr.debug(ve);
		}
	}
}
