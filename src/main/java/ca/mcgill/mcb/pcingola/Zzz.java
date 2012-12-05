package ca.mcgill.mcb.pcingola;

import ca.mcgill.mcb.pcingola.outputFormatter.VcfOutputFormatter;

/**
 * Simple test program
 * @author pcingola
 */
public class Zzz {

	public static void main(String[] args) {
		VcfOutputFormatter vof = new VcfOutputFormatter(null);
		String testIn[] = { "Hi ", "Hi how;", "Hi how;are|", "Hi how;are|you,", "Hi how;are|you,doing=", "Hi how;are|you,doing=today." };
		String testOut[] = { "Hi ", "Hi_how_", "Hi_how_are_", "Hi_how_are_you_", "Hi_how_are_you_doing_", "Hi_how_are_you_doing_today." };
		for (int i = 0; i < testIn.length; i++) {
			System.out.println("'" + testIn[i] + "'\t'" + vof.vcfInfoSafeString(testIn[i]) + "'\t'" + testOut[i] + "'");
		}

	}
}
