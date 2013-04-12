package ca.mcgill.mcb.pcingola;

/**
 * Tast class
 * @author pablocingolani
 */
public class Zzz {
	/**
	 * Main
	 * @param args
	 */
	public static void main(String[] args) {
		String str = "adsf asdf ; asdfasdf , asdfasfd = 5678578";
		String str2 = str.trim().replaceAll("(,|;|=| |\t)+", "_");
		System.out.println("STR2 = '" + str2 + "'");
	}

}
