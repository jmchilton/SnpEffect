package ca.mcgill.mcb.pcingola.snpEffect.commandLine;

import java.util.ArrayList;
import java.util.HashMap;

import ca.mcgill.mcb.pcingola.Config2DownloadTable;
import ca.mcgill.mcb.pcingola.Pcingola;
import ca.mcgill.mcb.pcingola.logStatsServer.LogStats;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * Command line program
 * 
 * @author pcingola
 */
public class SnpEff implements CommandLine {

	/**
	 *  Available gene database formats
	 */
	public enum GeneDatabaseFormat {
		BIOMART, GFF3, GFF2, GTF22, REFSEQ, GENBANK, EMBL
	}

	/**
	 *  Available input formats
	 */
	public enum InputFormat {
		TXT, PILEUP, VCF, BED
	}

	/**
	 *  Available output formats
	 */
	public enum OutputFormat {
		TXT, VCF, BED, BEDANN
	}

	public static final int COMMAND_LINE_WIDTH = 40;

	public static final String BUILD = "2012-04-20";
	public static final String VERSION_MAJOR = "2.1";
	public static final String REVISION = "b";
	public static final String VERSION_SHORT = VERSION_MAJOR + REVISION;
	public static final String VERSION = VERSION_SHORT + " (build " + BUILD + "), by " + Pcingola.BY;

	public static final String DEFAULT_SUMMARY_FILE = "snpEff_summary.html";
	public static final String DEFAULT_SUMMARY_GENES_FILE = "snpEff_genes.txt";

	// Parameters for LOG thread (a thread that logs information to a server)
	public static final int LOG_THREAD_WAIT_TIME = 1000; // 1 Second
	public static final int LOG_THREAD_WAIT_TIME_REPEAT = 3;

	protected String command = "";
	protected String[] args; // Arguments used to invoke this command
	protected String[] shiftArgs;
	protected boolean verbose = false; // Be verbose
	protected boolean quiet = false; // Be quiet
	protected boolean log = true; // Log to server (statistics)
	protected boolean multiThreaded = false; // Use multiple threads
	protected int numWorkers = Gpr.NUM_CORES; // Max number of threads (if multi-threaded version is available)
	protected int inOffset = 1; // By default positions are 1-based
	protected int outOffset = 1;
	protected String configFile; // Config file
	protected String genomeVer; // Genome version
	protected Config config; // Configuration

	/**
	 * Main
	 */
	public static void main(String[] args) {
		// Parse
		SnpEff snpEff = new SnpEff();
		snpEff.parseArgs(args);

		// Run
		boolean ok = snpEff.run();
		int retCode = ok ? 0 : -1;

		// Exit
		System.exit(retCode);
	}

	public SnpEff() {
		genomeVer = ""; // Genome version
		configFile = Config.DEFAULT_CONFIG_FILE; // Config file
	}

	/**
	 * 	Command line argument list (try to fit it into COMMAND_LINE_WIDTH)
	 * 
	 * @param splitLines
	 * @return
	 */
	String commandLineStr(boolean splitLines) {
		StringBuilder argsList = new StringBuilder();
		argsList.append("SnpEff " + command + " ");
		int size = argsList.length();

		for (String arg : args) {
			argsList.append(arg);
			size += arg.length();
			if (splitLines && (size > COMMAND_LINE_WIDTH)) {
				argsList.append(" \n");
				size = 0;
			} else {
				argsList.append(" ");
				size++;
			}
		}

		return argsList.toString();
	}

	/**
	 * Show an error (if not 'quiet' mode)
	 * @param message
	 */
	public void error(Throwable e, String message) {
		if (verbose && (e != null)) e.printStackTrace();
		if (!quiet) System.err.println(message);
	}

	/**
	 * Show an error message and exit
	 * @param message
	 */
	public void fatalError(String message) {
		System.err.println(message);
		System.exit(-1);
	}

	/**
	 * Parse command line arguments
	 * @param args
	 */
	@Override
	public void parseArgs(String[] args) {
		this.args = args;
		if (args.length <= 0) usage("Missing command");

		if (args[0].equalsIgnoreCase("build") //
				|| args[0].equalsIgnoreCase("dump") //
				|| args[0].equalsIgnoreCase("cds") //
				|| args[0].equalsIgnoreCase("eff") //
				|| args[0].equalsIgnoreCase("download") //
				|| args[0].equalsIgnoreCase("protein") //
				|| args[0].equalsIgnoreCase("closestExon") //
				|| args[0].equalsIgnoreCase("test") //
				|| args[0].equalsIgnoreCase("cfg2table") //
		) {
			command = args[0].toLowerCase();

			// Copy all args except initial 'command'
			ArrayList<String> argsList = new ArrayList<String>();
			for (int i = 1; i < args.length; i++) {
				if (args[i].equals("-noLog")) log = false; // This option is always available (to allow privacy in all commands)
				else argsList.add(args[i]);

				if (args[i].equals("-v")) verbose = true; // Make this option availabe here as well
			}
			shiftArgs = argsList.toArray(new String[0]);
		} else {
			command = "eff"; // Default command is 'eff'
			shiftArgs = args;
		}
	}

	/**
	 * Report stats to server
	 * @param ok
	 * @param errorMessage
	 */
	void report(boolean ok, String errorMessage, HashMap<String, String> reportValues) {
		if (!log) return; // Do nothing

		//---
		// Create logStats & add data 
		//---
		LogStats logStats = new LogStats();

		// Add paramters
		for (int i = 0; i < args.length; i++)
			logStats.add("args_" + i, args[i]);

		// Add run status info
		logStats.add("Finished_OK", Boolean.toString(ok));
		if (!errorMessage.isEmpty()) logStats.add("Error", errorMessage);

		// Add other info
		for (String name : reportValues.keySet())
			logStats.add(name, reportValues.get(name));

		//---
		// Run thread
		//---
		logStats.start();

		// Finish up
		if (verbose) Timer.showStdErr("Finishing up");
		for (int i = 0; i < LOG_THREAD_WAIT_TIME_REPEAT; i++) {
			if (!logStats.isAlive()) break;

			try {
				Thread.sleep(LOG_THREAD_WAIT_TIME); // Sleep 1 sec
			} catch (InterruptedException e) {
				; // Nothing to do
			}
		}

		// Interrupt?
		if (logStats.isAlive() && !logStats.isInterrupted()) {
			// Some people freak out about this 'Interrupting thread' message
			// if( verbose ) Timer.showStdErr("Interrupting thread");
			logStats.interrupt();
		}
	}

	/**
	 * Additional values to be reported
	 * @return
	 */
	public HashMap<String, String> reportValues() {
		HashMap<String, String> reportValues = new HashMap<String, String>();

		// Extra statistics: What kind of systems do users run this program on?
		String properties[] = { "user.name", "os.name", "os.version", "os.arch" };
		for (String prop : properties) {
			try {
				reportValues.put(prop, System.getProperty(prop));
			} catch (Exception e) {
				; // Do nothing, just skip the values
			};
		}

		try {
			reportValues.put("num.cores", Gpr.NUM_CORES + "");
			reportValues.put("total.mem", Runtime.getRuntime().totalMemory() + "");
		} catch (Exception e) {
			; // Do nothing, just skip the values
		};

		return reportValues;
	}

	/**
	 * Run according to command line options
	 */
	@Override
	public boolean run() {
		boolean ok = false;
		SnpEff snpEff = null;

		if (command.equals("build")) {
			//---
			// Build database
			//---
			snpEff = new SnpEffCmdBuild();
			snpEff.parseArgs(shiftArgs);
		} else if (command.equals("dump")) {
			//---
			// Dump database
			//---
			snpEff = new SnpEffCmdDump();
			snpEff.parseArgs(shiftArgs);
		} else if (command.equals("download")) {
			//---
			// Download database
			//---
			snpEff = new SnpEffCmdDownload();
			snpEff.parseArgs(shiftArgs);
		} else if (command.equals("cds")) {
			//---
			// CDS test
			//---
			snpEff = new SnpEffCmdCds();
			snpEff.parseArgs(shiftArgs);
		} else if (command.equals("eff")) {
			//---
			// Align to reference genome
			//---
			snpEff = new SnpEffCmdEff();
			snpEff.parseArgs(shiftArgs);
		} else if (command.equals("protein")) {
			//---
			// Protein test
			//---
			snpEff = new SnpEffCmdProtein();
			snpEff.parseArgs(shiftArgs);
		} else if (command.equals("closestexon")) {
			//---
			// Find closest exon
			//---
			snpEff = new SnpEffCmdClosestExon();
			snpEff.parseArgs(shiftArgs);
		} else if (command.equals("cfg2table")) {
			// Create download table and galaxy list from config file
			snpEff = new Config2DownloadTable();
			snpEff.parseArgs(shiftArgs);
		} else if (command.equals("test")) {
			//---
			// Test command (only for testing weird stuff)
			//---
			snpEff = new SnpEffCmdTest();
			snpEff.parseArgs(shiftArgs);
		} else throw new RuntimeException("Unknown command '" + command + "'");

		//---
		// Run
		//---
		String err = "";
		try {
			ok = snpEff.run();
		} catch (Throwable t) {
			err = t.getMessage();
			t.printStackTrace();
		}

		report(ok, err, snpEff.reportValues()); // Report 

		return ok;
	}

	/**
	 * Show 'usage' message and exit with an error code '-1'
	 * @param message
	 */
	@Override
	public void usage(String message) {
		if (message != null) System.err.println("Error: " + message + "\n");
		System.err.println("snpEff version " + VERSION);
		System.err.println("Usage: snpEff [eff]        [options] genome_version snp_file");
		System.err.println("   or: snpEff download     [options] genome_version");
		System.err.println("   or: snpEff build        [options] genome_version");
		System.err.println("   or: snpEff dump         [options] genome_version");
		System.err.println("   or: snpEff cds          [options] genome_version");
		System.err.println("   or: snpEff protein      [options] genome_version");
		System.err.println("   or: snpEff closestExon  [options] genome_version");
		System.exit(-1);
	}
}
