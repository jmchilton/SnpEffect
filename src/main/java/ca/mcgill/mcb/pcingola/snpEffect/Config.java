package ca.mcgill.mcb.pcingola.snpEffect;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Properties;

import ca.mcgill.mcb.pcingola.codons.CodonTable;
import ca.mcgill.mcb.pcingola.codons.CodonTables;
import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.interval.Genome;
import ca.mcgill.mcb.pcingola.util.Gpr;

@SuppressWarnings("serial")
public class Config implements Serializable, Iterable<String> {

	public static final String DEFAULT_CONFIG_FILE = "snpEff.config";
	public static final String DEFAULT_DATA_DIR = "./data";
	public static String GENOMES_DIR = "genomes"; // Directory has one genomes information (FASTA files)

	public static final String GENOME_KEY = ".genome";
	public static final String REFERENCE_KEY = ".reference";
	public static final String CODON_KEY = "codon.";
	public static final String CODONTABLE_KEY = ".codonTable";

	public static boolean debug = false; // Debug mode?
	private static Config configInstance = null; // Config is some kind of singleton because we want to make it accessible from everywhere

	boolean treatAllAsProteinCoding = true; // Calculate effect only in coding genes. Default is true for testing and debugging reasons (command line default is 'false')
	boolean onlyRegulation = false; // Only use regulation features
	boolean errorOnMissingChromo = true; // Error if chromosome is missing
	boolean errorChromoHit = true; // Error if chromosome is not hit in a query
	String dataDir;
	Genome genome;
	HashMap<String, Genome> genomeByVersion;
	HashMap<String, String> referenceByVersion;
	HashMap<String, String> nameByVersion;
	SnpEffectPredictor snpEffectPredictor;
	String databaseRepository = "";

	public static Config get() {
		return configInstance;
	}

	/**
	 * Create a config (uses DEFAULT_CONFIG_FILE)
	 * @param genomeVersion
	 */
	public Config(String genomeVersion) {
		read(genomeVersion, DEFAULT_CONFIG_FILE); // Read config file and get a genome
		genome = genomeByVersion.get(genomeVersion); // Set a genome
		if (genome == null) throw new RuntimeException("No such genome '" + genomeVersion + "'");
		configInstance = this;
	}

	/**
	 * Create a configuration from 'configFileName'
	 * @param genomeVersion
	 * @param configFileName
	 */
	public Config(String genomeVersion, String configFileName) {
		read(genomeVersion, configFileName); // Read config file and get a genome
		genome = genomeByVersion.get(genomeVersion); // Set a genome
		if (genome == null) throw new RuntimeException("No such genome '" + genomeVersion + "'");
		configInstance = this;
	}

	/**
	 * Extract and create codon tables
	 * @param genomeVersion
	 * @param properties
	 */
	void createCodonTables(String genomeVersion, Properties properties) {
		//---
		// Read codon tables
		//---
		for (Object key : properties.keySet()) {
			if (key.toString().startsWith(CODON_KEY)) {
				String name = key.toString().substring(CODON_KEY.length());
				String table = properties.getProperty(key.toString());
				CodonTable codonTable = new CodonTable(name, table);
				CodonTables.getInstance().add(codonTable);
			}
		}

		//---
		// Assign codon tables for different genome+chromosome
		//---
		for (Object key : properties.keySet()) {
			String keyStr = key.toString();
			if (keyStr.endsWith(CODONTABLE_KEY) && keyStr.startsWith(genomeVersion)) {
				// We have to find genome name and chromosome name
				// E.g. parse "dm5.30.dmel_mitochondrion_genome.codonTable" => genomeVer="dm5.30" and chromo="dmel_mitochondrion_genome"
				String genomeVer = findGenome(keyStr);
				if (genomeVer == null) throw new RuntimeException("Cannot find appropriate genome for '" + keyStr + "'");

				// Everything between gneomeName and ".codonTable" is assumed to be chromosome name
				String chromo = keyStr.substring(genomeVer.length() + 1, keyStr.length() - CODONTABLE_KEY.length());

				// Ignore configuration for other genomes
				if (genomeVersion.equals(genomeVer)) {

					String codonTableName = properties.getProperty(key.toString());
					CodonTable codonTable = CodonTables.getInstance().getTable(codonTableName);

					// Sanity checks
					if (genomeByVersion.get(genomeVer) == null) throw new RuntimeException("Error parsing property '" + key + "'. No such genome '" + genomeVer + "'");
					if (!genomeByVersion.get(genomeVer).hasChromosome(chromo)) {
						// Create chromosome
						Genome genome = genomeByVersion.get(genomeVer);
						Chromosome chr = new Chromosome(genome, 0, 0, 1, chromo);
						genome.add(chr);
					}
					if (codonTable == null) throw new RuntimeException("Error parsing property '" + key + "'. No such codon table '" + codonTableName + "'");

					// Everything seems to be OK, go on
					CodonTables.getInstance().add(genomeByVersion.get(genomeVer), chromo, codonTable);
				}
			}
		}
	}

	/**
	 * Find the genome string in a 'codonTable' key string
	 * @param codonKeyStr
	 * @return
	 */
	String findGenome(String codonKeyStr) {
		for (String genVer : genomeByVersion.keySet()) {
			if (codonKeyStr.startsWith(genVer + ".")) return genVer;
		}
		return null;
	}

	/**
	 * Genes file path (no extension)
	 * @return
	 */
	public String getBaseFileNameGenes() {
		return dataDir + "/" + genome.getVersion() + "/genes";
	}

	/**
	 * Regulation file (GFF format)
	 * @return
	 */
	public String getBaseFileNameRegulation() {
		return getDirDataVersion() + "/regulation";
	}

	public String getDatabaseRepository() {
		return databaseRepository;
	}

	/**
	 * Main data directory 
	 * @return
	 */
	public String getDirData() {
		return dataDir;
	}

	/**
	 * Data dir for a specific genome version (i.e. where the database is) 
	 * @return
	 */
	public String getDirDataVersion() {
		return dataDir + "/" + genome.getVersion();
	}

	/**
	 * Directory where regulation 'BED' files are
	 * @return
	 */
	public String getDirRegulationBed() {
		return getDirDataVersion() + "/regulation.bed/";
	}

	public String getFileNameCds() {
		return getDirDataVersion() + "/cds.fa";
	}

	/**
	 * Filenames for reference sequence (fasta files)
	 * @return
	 */
	public List<String> getFileNameGenomeFasta() {
		ArrayList<String> files = new ArrayList<String>();
		files.add(getDirData() + "/genomes/" + genome.getVersion() + ".fa");
		files.add(getDirData() + "/" + genome.getVersion() + "/sequences.fa");
		return files;
	}

	public String getFileNameProteins() {
		return getDirDataVersion() + "/protein.fa";
	}

	public Genome getGenome() {
		return genome;
	}

	public String getName(String genomeVersion) {
		return nameByVersion.get(genomeVersion);
	}

	public String getReference(String genomeVersion) {
		return referenceByVersion.get(genomeVersion);
	}

	public SnpEffectPredictor getSnpEffectPredictor() {
		return snpEffectPredictor;
	}

	public boolean isErrorChromoHit() {
		return errorChromoHit;
	}

	public boolean isErrorOnMissingChromo() {
		return errorOnMissingChromo;
	}

	public boolean isOnlyRegulation() {
		return onlyRegulation;
	}

	public boolean isTreatAllAsProteinCoding() {
		return treatAllAsProteinCoding;
	}

	@Override
	public Iterator<String> iterator() {
		return nameByVersion.keySet().iterator();
	}

	/**
	 * Load a snpEff predictor
	 * WARNING: 'genome' object get replaced upon loading a snpEffectPredictor (this is a dangerous side effect)
	 */
	public SnpEffectPredictor loadSnpEffectPredictor() {
		snpEffectPredictor = SnpEffectPredictor.load(this);
		genome = snpEffectPredictor.genome; // WARNING: 'genome' object get replaced upon loading a snpEffectPredictor (this might have dangerous side effects)
		return snpEffectPredictor;
	}

	/**
	 * Read configuration file and create all 'genomes' 
	 * @return
	 */
	private void read(String genomeVersion, String configFileName) {
		//---
		// Read properties file
		//---
		Properties properties = new Properties();
		try {
			properties.load(new FileReader(new File(configFileName)));
		} catch (FileNotFoundException e) {
			throw new RuntimeException("Cannot find config file '" + configFileName + "'");
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		// Config directory
		String cfgDir = "";
		try {
			File configDir = new File(configFileName).getAbsoluteFile().getParentFile();
			cfgDir = configDir.getCanonicalPath();
		} catch (IOException e1) {
		}

		//---
		// Set attributes
		//---
		dataDir = properties.getProperty("data_dir", DEFAULT_DATA_DIR);
		if (dataDir.startsWith("~")) dataDir = Gpr.HOME + "/" + dataDir.substring(1); // Relative to 'home' dir?
		else if (!dataDir.startsWith("/")) dataDir = cfgDir + "/" + dataDir; // Not an absolute path?
		if (dataDir.endsWith("/")) dataDir = dataDir.substring(0, dataDir.length() - 1); // make sure path doesn't end with '/' (some OS can have problems with "//" in paths)

		databaseRepository = properties.getProperty("database_repository", "");

		//---
		// Find all genomes in this config file
		//---
		genomeByVersion = new HashMap<String, Genome>();
		referenceByVersion = new HashMap<String, String>();
		nameByVersion = new HashMap<String, String>();
		for (Object k : properties.keySet()) {
			String key = k.toString();
			if (key.endsWith(GENOME_KEY)) {
				String genVer = key.substring(0, key.length() - GENOME_KEY.length());

				// Add full namne
				String name = properties.getProperty(genVer + GENOME_KEY);
				nameByVersion.put(genVer, name);

				// Add reference
				String ref = properties.getProperty(genVer + REFERENCE_KEY);
				referenceByVersion.put(genVer, ref);
			}
		}

		// Read config file for genome version (if any)
		readGenomeConfig(genomeVersion, properties);

		// Codon tables
		createCodonTables(genomeVersion, properties);
	}

	/**
	 * Read a config file for a given genome version (dataDir/genVer/snpEff.config)
	 * Add all properties found to 'properties'
	 * 
	 * @param genVer
	 * @param properties
	 */
	void readGenomeConfig(String genVer, Properties properties) {
		String genomePropsFileName = dataDir + "/" + genVer + "/snpEff.config";
		try {
			// Read properties file "data_dir/genomeVersion/snpEff.conf"
			// If the file exists, read all properties and add them to the main 'properties' 
			Properties genomeProps = new Properties();
			genomeProps.load(new FileReader(new File(genomePropsFileName)));

			// Copy genomeProperties to 'properties'
			for (Object propKey : genomeProps.keySet()) {
				String pk = propKey.toString();
				String propVal = genomeProps.getProperty(pk);
				if (properties.getProperty(pk) == null) {
					properties.setProperty(pk, propVal);
				} else System.err.println("Ignoring property '" + pk + "' in file '" + genomePropsFileName + "'");
			}
		} catch (Exception e) {
			if (debug) System.err.println("File '" + genomePropsFileName + "' not found"); // File does not exists? => OK
		}

		Genome genome = new Genome(genVer, properties);
		genomeByVersion.put(genVer, genome);
	}

	public void setErrorChromoHit(boolean errorChromoHit) {
		this.errorChromoHit = errorChromoHit;
	}

	public void setErrorOnMissingChromo(boolean errorOnMissingChromo) {
		this.errorOnMissingChromo = errorOnMissingChromo;
	}

	public void setOnlyRegulation(boolean onlyRegulation) {
		this.onlyRegulation = onlyRegulation;
	}

	public void setSnpEffectPredictor(SnpEffectPredictor snpEffectPredictor) {
		this.snpEffectPredictor = snpEffectPredictor;
	}

	public void setTreatAllAsProteinCoding(boolean treatAllAsProteinCoding) {
		this.treatAllAsProteinCoding = treatAllAsProteinCoding;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();

		for (String genVer : this) {
			String name = nameByVersion.get(genVer).replace('_', ' ');
			String ref = referenceByVersion.get(genVer);

			sb.append("\t" + genVer);
			sb.append("\t" + name);
			if (ref != null) sb.append("\t" + ref);
			sb.append("\n");
		}
		return sb.toString();
	}
}
