package ca.mcgill.mcb.pcingola.snpEffect.commandLine;

import java.util.ArrayList;
import java.util.HashMap;

import ca.mcgill.mcb.pcingola.Pcingola;
import ca.mcgill.mcb.pcingola.interval.Custom;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Markers;
import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.logStatsServer.LogStats;
import ca.mcgill.mcb.pcingola.logStatsServer.VersionCheck;
import ca.mcgill.mcb.pcingola.nextProt.SnpEffCmdBuildNextProt;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.spliceSites.SpliceAnalysis;
import ca.mcgill.mcb.pcingola.util.Gpr;

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
		BIOMART, GFF3, GFF2, GTF22, REFSEQ, KNOWN_GENES, GENBANK, EMBL
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
		TXT, VCF, BED, BEDANN, GATK
	}

	public static final int COMMAND_LINE_WIDTH = 40;

	public static final String SOFTWARE_NAME = "SnpEff";
	public static final String REVISION = "";
	public static final String BUILD = "2013-05-14";
	public static final String VERSION_MAJOR = "3.2";
	public static final String VERSION_SHORT = VERSION_MAJOR + REVISION;
	public static final String VERSION_NO_NAME = VERSION_SHORT + " (build " + BUILD + "), by " + Pcingola.BY;
	public static final String VERSION = SOFTWARE_NAME + " " + VERSION_NO_NAME;
	public static final String DEFAULT_SUMMARY_FILE = "snpEff_summary.html";
	public static final String DEFAULT_SUMMARY_GENES_FILE = "snpEff_genes.txt";

	protected String command = "";
	protected String[] args; // Arguments used to invoke this command
	protected String[] shiftArgs;
	protected boolean help; // Show command help and exit
	protected boolean verbose; // Be verbose
	protected boolean debug; // Debug mode
	protected boolean quiet; // Be quiet
	protected boolean log; // Log to server (statistics)
	protected boolean multiThreaded; // Use multiple threads
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
		System.exit(ok ? 0 : -1);
	}

	public SnpEff() {
		genomeVer = ""; // Genome version
		configFile = Config.DEFAULT_CONFIG_FILE; // Config file
		verbose = false; // Be verbose
		debug = false; // Debug mode
		quiet = false; // Be quiet
		log = true; // Log to server (statistics)
		multiThreaded = false; // Use multiple threads
	}

	/**
	 * Check if there is a new version of the program
	 */
	void checkNewVersion(Config config) {
		// Download command checks for versions, no need to do it twice
		if ((config != null) && !command.equalsIgnoreCase("download")) {
			// Check if a new version is available
			VersionCheck versionCheck = VersionCheck.version(SnpEff.SOFTWARE_NAME, SnpEff.VERSION_SHORT, config.getVersionsUrl(), verbose);
			if (!quiet && versionCheck.isNewVersion()) {
				System.err.println("\n\nNEW VERSION!\n\tThere is a new " + this.getClass().getSimpleName() + " version available: " //
						+ "\n\t\tVersion      : " + versionCheck.getLatestVersion() // 
						+ "\n\t\tRelease date : " + versionCheck.getLatestReleaseDate() //
						+ "\n\t\tDownload URL : " + versionCheck.getLatestUrl() //
						+ "\n" //
				);
			}
		}
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
	 * Is this a command line option (e.g. "-tfam" is a command line option, but "-" means STDIN)
	 * @param arg
	 * @return
	 */
	protected boolean isOpt(String arg) {
		return arg.startsWith("-") && (arg.length() > 1);
	}

	/**
	 * Parse command line arguments
	 * @param args
	 */
	@Override
	public void parseArgs(String[] args) {
		this.args = args;
		if (args.length <= 0) usage("Missing command");

		int argNum = 0;

		//---
		// Parse command
		//---
		if (args[0].equalsIgnoreCase("build") //
				|| args[0].equalsIgnoreCase("buildNextProt") //
				|| args[0].equalsIgnoreCase("dump") //
				|| args[0].equalsIgnoreCase("cds") //
				|| args[0].equalsIgnoreCase("eff") //
				|| args[0].equalsIgnoreCase("download") //
				|| args[0].equalsIgnoreCase("protein") //
				|| args[0].equalsIgnoreCase("closest") //
				|| args[0].equalsIgnoreCase("test") //
				|| args[0].equalsIgnoreCase("databases") //
				|| args[0].equalsIgnoreCase("spliceAnalysis") //
				|| args[0].equalsIgnoreCase("count") //
				|| args[0].equalsIgnoreCase("genes2bed") //
				|| args[0].equalsIgnoreCase("len") //
		) {
			command = args[argNum++].toLowerCase();
		} else {
			command = "eff"; // Default command is 'eff'
		}

		//---
		// Copy all args except initial 'command'
		//---
		ArrayList<String> argsList = new ArrayList<String>();
		for (int i = argNum; i < args.length; i++) {
			// These options are available for allow all commands
			if (args[i].equalsIgnoreCase("-noLog")) log = false;
			else if (args[i].equals("-h") || args[i].equalsIgnoreCase("-help")) help = true;
			else if (args[i].equals("-v") || args[i].equalsIgnoreCase("-verbose")) {
				verbose = true;
				quiet = false;
			} else if (args[i].equals("-q") || args[i].equalsIgnoreCase("-quiet")) {
				quiet = true;
				verbose = false;
			} else if (args[i].equals("-d") || args[i].equalsIgnoreCase("-debug")) debug = verbose = true;
			else if ((args[i].equals("-c") || args[i].equalsIgnoreCase("-config"))) {
				if ((i + 1) < args.length) configFile = args[++i];
				else usage("Option '-c' without config file argument");
			} else argsList.add(args[i]);
		}

		shiftArgs = argsList.toArray(new String[0]);
	}

	/**
	 * Read markers file
	 * 
	 * Supported formats: BED, TXT, BigBed
	 * 
	 * @param fileName
	 * @return
	 */
	protected Markers readMarkers(String fileName) {
		Markers markersSeqChange = Markers.readMarkers(fileName);
		String label = Gpr.removeExt(Gpr.baseName(fileName));

		// Convert 'SeqChange' markers to 'Custom' markers
		Markers markers = new Markers();
		for (Marker m : markersSeqChange) {
			Custom custom = new Custom(m.getParent(), m.getStart(), m.getEnd(), m.getStrand(), m.getId(), label);
			custom.setScore(((SeqChange) m).getScore());
			markers.add(custom);
		}

		// Number added
		return markers;
	}

	/**
	 * Additional values to be reported
	 * @return
	 */
	public HashMap<String, String> reportValues() {
		HashMap<String, String> reportValues = new HashMap<String, String>();
		return reportValues;
	}

	/**
	 * Run according to command line options
	 */
	@Override
	public boolean run() {
		boolean ok = false;
		SnpEff snpEff = null;

		// All commands are lowercase
		command = command.toLowerCase();
		if (command.equalsIgnoreCase("build")) snpEff = new SnpEffCmdBuild();
		else if (command.equalsIgnoreCase("buildNextProt")) snpEff = new SnpEffCmdBuildNextProt();
		else if (command.equalsIgnoreCase("dump")) snpEff = new SnpEffCmdDump();
		else if (command.equalsIgnoreCase("download")) snpEff = new SnpEffCmdDownload();
		else if (command.equalsIgnoreCase("cds")) snpEff = new SnpEffCmdCds();
		else if (command.equalsIgnoreCase("eff")) snpEff = new SnpEffCmdEff();
		else if (command.equalsIgnoreCase("protein")) snpEff = new SnpEffCmdProtein();
		else if (command.equalsIgnoreCase("closest")) snpEff = new SnpEffCmdClosest();
		else if (command.equalsIgnoreCase("databases")) snpEff = new SnpEffCmdDatabases();
		else if (command.equalsIgnoreCase("genes2bed")) snpEff = new SnpEffCmdGenes2Bed();
		else if (command.equalsIgnoreCase("spliceanalysis")) snpEff = new SpliceAnalysis();
		else if (command.equalsIgnoreCase("count")) snpEff = new SnpEffCmdCount();
		else if (command.equalsIgnoreCase("len")) snpEff = new SnpEffCmdLen();
		else if (command.equalsIgnoreCase("test")) snpEff = new SnpEffCmdTest();
		else throw new RuntimeException("Unknown command '" + command + "'");

		//---
		// Run
		//---
		String err = "";
		try {
			snpEff.verbose = verbose;
			snpEff.help = help;
			snpEff.debug = debug;
			snpEff.quiet = quiet;
			snpEff.configFile = configFile;

			if (help) snpEff.usage(null); // Show help message and exit
			else snpEff.parseArgs(shiftArgs);
			ok = snpEff.run();
		} catch (Throwable t) {
			err = t.getMessage();
			t.printStackTrace();
		}

		// Report to server (usage statistics) 
		if (log) {
			// Log to server
			LogStats.report(SOFTWARE_NAME, VERSION_SHORT, VERSION, ok, verbose, args, err, snpEff.reportValues());

			// Check for new version (use config file from command, since this one doesn't load a config file)
			checkNewVersion(snpEff.config);
		}

		return ok;
	}

	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}

	/**
	 * Show 'usage' message and exit with an error code '-1'
	 * @param message
	 */
	@Override
	public void usage(String message) {
		if (message != null) System.err.println("Error: " + message + "\n");
		System.err.println("snpEff version " + VERSION);
		System.err.println("Usage: snpEff [command] [options] [files]");
		System.err.println("\nAvailable commands: ");
		System.err.println("   [eff]           : Calculate effect of variants. Default: eff (no command or 'eff').");
		System.err.println("   build           : Build a SnpEff database.");
		System.err.println("   buildNextProt   : Build a SnpEff for NextProt (using NextProt's XML files).");
		System.err.println("   cds             : Compare CDS sequences calculated form a SnpEff database to the one in a FASTA file. Used for checking databases correctness.");
		System.err.println("   closest         : Annotate the closest genomic region.");
		System.err.println("   count           : Count how many intervals (from a BAM, BED or VCF file) overlap with each genomic interval.");
		System.err.println("   databases       : Show currently available databases (from local config file).");
		System.err.println("   download        : Download a SnpEff database.");
		System.err.println("   dump            : Dump to STDOUT a SnpEff database (mostly used for debugging).");
		System.err.println("   genes2bed       : Create a bed file from a genes list.");
		System.err.println("   len             : Calculate total genomic length for each marker type.");
		System.err.println("   protein         : Compare protein sequences calculated form a SnpEff database to the one in a FASTA file. Used for checking databases correctness.");
		System.err.println("   spliceAnalysis  : Perform an analysis of splice sites. Experimental feature.");
		System.err.println("\nRun 'java -jar snpEff.jar command' for help on each specific command");
		System.exit(-1);
	}
}
