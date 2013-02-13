package ca.mcgill.mcb.pcingola.snpEffect.commandLine;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import ca.mcgill.mcb.pcingola.snpEffect.Config;

/**
 * Show all databases configures in snpEff.config
 * 
 * Create an HTML 'download' table based on the config file
 * Also creates a list of genome for Galaxy menu
 * 
 * @author pablocingolani
 */
public class SnpEffCmdDatabases extends SnpEff {

	public static final String DARK_ROW = "bgcolor=#CCCCCC";
	public static final String LIGHT_ROW = "bgcolor=#EEEEEE";

	public static final String HTTP_PROTOCOL = "http://";
	public static final String FTP_PROTOCOL = "ftp://";

	boolean galaxy = false;
	boolean html = false;
	Config config;
	HashMap<String, String> nameByGenomeVer;
	ArrayList<String> namesSorted;
	ArrayList<String> genVerSorted;

	public SnpEffCmdDatabases() {
		// Read config (it doesn't matter which genome)
		config = new Config("hg19");

		// Get all genome names and sort them
		nameByGenomeVer = new HashMap<String, String>();
		for (String genVer : config)
			nameByGenomeVer.put(genVer, config.getName(genVer));

		namesSorted = new ArrayList<String>();
		namesSorted.addAll(nameByGenomeVer.values());
		Collections.sort(namesSorted);

		// Sort genome versions
		genVerSorted = new ArrayList<String>();
		for (String genVer : config)
			genVerSorted.add(genVer);
		Collections.sort(genVerSorted);
	}

	/**
	 * Galaxy config genome list
	 */
	void galaxyConfig() {
		System.out.println("\t<param name=\"genomeVersion\" type=\"select\" label=\"Genome\">");

		for (String name : namesSorted) {
			for (String genVer : genVerSorted) {
				String n = config.getName(genVer);

				// In this group?
				if (name.equals(n)) {
					System.out.println("\t\t<option value=\"" + genVer + "\">" + name.replace('_', ' ') + " : " + genVer + "</option>");
				}
			}
		}

		System.out.println("\t</param>");
	}

	/**
	 * Create html table
	 */
	void htmlTable() {
		// Create an HTML table
		boolean dark = false;
		String bg = "";

		System.out.println("\t<table> <tr " + DARK_ROW + "> <td> <b> Genome </b> </td>  <td> <b> Version </b> </td>  <td> <b> Reference </b> </td> </tr>");
		for (String name : namesSorted) {

			// Color
			if (dark) bg = DARK_ROW;
			else bg = LIGHT_ROW;
			dark = !dark;

			boolean showName = true;
			for (String genVer : genVerSorted) {
				String n = config.getName(genVer);
				// In this group?
				if (name.equals(n)) {
					System.out.println("\t\t<tr " + bg + ">");

					// Show name
					String name2show = showName ? name.replace('_', ' ') : "&nbsp;";
					System.out.println("\t\t\t<td> " + name2show + " </td>");
					showName = false;

					// Download link
					String url = "http://sourceforge.net/projects/snpeff/files/databases/v" + SnpEff.VERSION_MAJOR + "/snpEff_v" + SnpEff.VERSION_MAJOR + "_" + genVer + ".zip";
					System.out.println("\t\t\t<td> <a class=\"body\" href=\"" + url + "\"> " + genVer + " </a> </td>");

					// Reference
					String ref = config.getReference(genVer);
					String link = "";
					if (ref != null) {
						if (ref.indexOf(',') > 0) ref = ref.substring(0, ref.indexOf(',')); // Many references? Use the first one
						link = ref;

						// Remove everything after slash
						int idx = ref.indexOf('/', HTTP_PROTOCOL.length());
						if (idx > 0) ref = ref.substring(0, idx);

						// Remove protocol
						if (ref.startsWith(HTTP_PROTOCOL)) ref = ref.substring(HTTP_PROTOCOL.length());
						if (ref.startsWith(FTP_PROTOCOL)) ref = ref.substring(FTP_PROTOCOL.length());

					} else ref = "";
					System.out.println("\t\t\t<td> <a class=\"body\" href=\"" + link + "\">" + ref + "</a> </td>");

					System.out.println("\t\t</tr>");
				}
			}
		}
		System.out.println("\t</table>");
	}

	@Override
	public void parseArgs(String[] args) {
		if (args.length > 1) usage(null);

		if (args.length == 1) {
			galaxy = args[0].equals("galaxy");
			html = args[0].equals("html");
		}
	}

	@Override
	public boolean run() {
		if (galaxy) galaxyConfig();
		else if (html) htmlTable();
		else txtTable();

		return true;
	}

	/**
	 * Create TXT table
	 */
	void txtTable() {
		System.out.println(String.format("%-60s\t%-60s\t%s", "Genome", "Organism", "Database download link"));
		System.out.println(String.format("%-60s\t%-60s\t%s", "------", "--------", "----------------------"));

		for (String genomeVer : genVerSorted) {
			String name = nameByGenomeVer.get(genomeVer);

			// Download link
			String url = "http://sourceforge.net/projects/snpeff/files/databases/v" + SnpEff.VERSION_MAJOR + "/snpEff_v" + SnpEff.VERSION_MAJOR + "_" + name + ".zip";

			// Show
			System.out.println(String.format("%-60s\t%-60s\t%s", genomeVer, name, url));
		}
	}

	@Override
	public void usage(String message) {
		if (message != null) System.err.println("Error: " + message + "\n");
		System.err.println("Usage: snpEff databases [galaxy]");
		System.err.println("\nOptions");
		System.err.println("\n\tgalaxy\t: Show databases in a galaxy menu format.");
		System.exit(-1);
	}
}