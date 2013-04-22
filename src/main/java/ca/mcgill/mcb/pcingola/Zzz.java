package ca.mcgill.mcb.pcingola;

public class Zzz {

	public static void main(String[] args) {
		String str = "asdfasdf(asdfadf]a  asdf asdf asdfasdf;fasdf,asdfad|asdfadf}".replaceAll("(\\s|\\(|\\)|\\[|\\]|;|,|\\|)+", "_");
		System.out.println(str);
	}
}
