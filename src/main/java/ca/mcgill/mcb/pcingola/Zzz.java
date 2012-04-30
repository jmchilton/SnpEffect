package ca.mcgill.mcb.pcingola;

public class Zzz {

	public static void main(String[] args) {
		String s = "Hello\n\tworld!";
		String ss = s.replaceAll("\\n", "\\\\n");
		ss = ss.replaceAll("\\t", "\\\\t");
		ss = ss.replaceAll("\\r", "");
		System.out.println(s);
		System.out.println("--------------------------------------");
		System.out.println(ss);
	}
}
