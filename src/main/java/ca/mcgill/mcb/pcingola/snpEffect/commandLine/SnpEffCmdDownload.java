package ca.mcgill.mcb.pcingola.snpEffect.commandLine;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.util.Date;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * Command line program: Build database
 * 
 * @author pcingola
 */
public class SnpEffCmdDownload extends SnpEff {

	public static boolean debug = false;

	String version = SnpEff.VERSION_MAJOR;

	private static int BUFFER_SIZE = 10480;

	public SnpEffCmdDownload() {
		super();
	}

	/**
	 * File name from URL (i.e. anything after the last '/')
	 * @param url
	 * @return
	 */
	String baseName(String url) {
		String f[] = url.toString().split("/");
		return f[f.length - 1];
	}

	/**
	 * Build the URL for getting the database file
	 * 
	 * Format  : DatabaseRepository / v VERSION / snpEff_v VERSION _ genomeVersion .zip
	 * Example : http://downloads.sourceforge.net/project/snpeff/databases/v2_0_3/snpEff_v2_0_3_EF3.64.zip
	 * 
	 * @param genomeVer
	 * @return
	 */
	private URL buildUrl() {
		try {
			// Replace '.' by '_' 
			version = version.replace('.', '_');

			String urlRoot = config.getDatabaseRepository();

			StringBuilder urlsb = new StringBuilder();
			urlsb.append(urlRoot);
			if (urlsb.charAt(urlRoot.length() - 1) != '/') urlsb.append("/");
			urlsb.append("v" + version + "/snpEff_v" + version + "_" + genomeVer + ".zip");

			return new URL(urlsb.toString());
		} catch (MalformedURLException e) {
			return null;
		}
	}

	/**
	 * Download a file
	 * @param url
	 * @return
	 */
	boolean download(URL url, String localFile) {
		boolean res = false;
		try {
			Timer.show("Connecting to " + url);
			URLConnection connection = url.openConnection();

			for (boolean followRedirect = true; followRedirect;) {
				HttpURLConnection httpConnection = (HttpURLConnection) connection;
				int code = httpConnection.getResponseCode();

				if (code == 200) {
					followRedirect = false; // We are done
				} else if (code == 302) {
					String newUrl = connection.getHeaderField("Location");
					Timer.show("Following redirect: " + newUrl);
					url = new URL(newUrl);
					connection = url.openConnection();
				} else if (code == 404) {
					throw new RuntimeException("File not found on the server. Make sure the database name is correct.");
				} else throw new RuntimeException("Error code from server: " + code);
			}

			// Copy resource to local file, use remote file if no local file name specified
			InputStream is = url.openStream();

			// Print info about resource
			Date date = new Date(connection.getLastModified());
			Timer.show("Copying file (type: " + connection.getContentType() + ", modified on: " + date + ")");

			// Open local file 
			Timer.show("Local file name: '" + localFile + "'");
			FileOutputStream os = null;
			os = new FileOutputStream(localFile);

			// Copy to file
			int count = 0, total = 0, lastShown = 0;
			byte data[] = new byte[BUFFER_SIZE];
			while ((count = is.read(data, 0, BUFFER_SIZE)) != -1) {
				os.write(data, 0, count);
				total += count;

				// Show every MB
				if ((total - lastShown) > (1024 * 1024)) {
					Timer.show("Downloaded " + total + " bytes");
					lastShown = total;
				}
			}

			// Close streams
			is.close();
			os.close();
			Timer.show("Donwload finished. Total " + total + " bytes.");

			res = true;
		} catch (Exception e) {
			throw new RuntimeException(e);
		}

		return res;
	}

	/**
	 * Parse command line arguments
	 * @param args
	 */
	@Override
	public void parseArgs(String[] args) {
		this.args = args;
		for (int i = 0; i < args.length; i++) {

			// Argument starts with '-'?
			if (args[i].startsWith("-")) {
				if ((args[i].equals("-c") || args[i].equalsIgnoreCase("-config"))) {
					if ((i + 1) < args.length) configFile = args[++i];
					else usage("Option '-c' without config file argument");
				} else if (args[i].equals("-v") || args[i].equalsIgnoreCase("-verbose")) {
					verbose = true;
					quiet = false;
				} else if (args[i].equals("-h") || args[i].equalsIgnoreCase("-help")) {
					usage(null);
					System.exit(0);
				} else usage("Unknow option '" + args[i] + "'");
			} else if (genomeVer.length() <= 0) genomeVer = args[i];
			else usage("Unknow parameter '" + args[i] + "'");
		}

		// Check: Do we have all required parameters?
		if (genomeVer.isEmpty()) usage("Missing genomer_version parameter");
	}

	/**
	 * Download database from server
	 */
	@Override
	public boolean run() {
		config = new Config(genomeVer, configFile);

		Timer.show("Downloading database for '" + genomeVer + "'");

		URL url = buildUrl();
		String localFile = baseName(url.toString());

		// Download and unzip
		if (download(url, localFile)) {
			if (unzip(localFile)) Timer.show("Unzip: OK");
		}

		Timer.show("Done");
		return true;
	}

	/**
	 * Unzip all files
	 * @return
	 */
	boolean unzip(String zipFile) {
		try {
			FileInputStream fis = new FileInputStream(zipFile);
			ZipInputStream zin = new ZipInputStream(new BufferedInputStream(fis));

			ZipEntry entry;
			while ((entry = zin.getNextEntry()) != null) {

				byte data[] = new byte[BUFFER_SIZE];
				Timer.show("Extracting file '" + entry.getName() + "'");

				//---
				// Move to 'data' dir
				//---
				String entryName = entry.getName();
				String entryPath[] = entryName.split("/"); // Entry name should be something like 'data/genomeVer/file';
				String dataName = entryPath[entryPath.length - 2] + "/" + entryPath[entryPath.length - 1]; // remove the 'data/' part
				entryName = config.getDirData() + "/" + dataName; // Ad local 'data' dir
				Timer.show("Local file name: '" + entryName + "'");

				// Create local dir
				String dir = Gpr.dirName(entryName);
				Timer.show("Creating local directory: '" + dir + "'");
				new File(dir).mkdirs();

				//---
				// Extract data
				//---
				FileOutputStream fos = new FileOutputStream(entryName);
				BufferedOutputStream dest = new BufferedOutputStream(fos, BUFFER_SIZE);

				int count = 0;
				while ((count = zin.read(data, 0, BUFFER_SIZE)) != -1)
					dest.write(data, 0, count);

				dest.flush();
				dest.close();
			}
			zin.close();

		} catch (Exception e) {
			throw new RuntimeException(e);
		}

		return true;
	}

	/**
	 * Show 'usage;' message and exit with an error code '-1'
	 * @param message
	 */
	@Override
	public void usage(String message) {
		if (message != null) System.err.println("Error: " + message + "\n");
		System.err.println("snpEff version " + VERSION);
		System.err.println("Usage: snpEff download [options] genome_version");
		System.err.println("\nGeneric options:");
		System.err.println("\t-c , -config            : Specify config file");
		System.err.println("\t-h , -help              : Show this help and exit");
		System.err.println("\t-v , -verbose           : Verbose mode");
		System.err.println("\t-noLog                  : Do not report usage statistics to server");
		System.exit(-1);
	}
}
