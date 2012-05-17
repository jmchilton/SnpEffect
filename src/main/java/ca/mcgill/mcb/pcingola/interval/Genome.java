package ca.mcgill.mcb.pcingola.interval;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Properties;

import ca.mcgill.mcb.pcingola.fileIterator.FastaFileIterator;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * 
 * This is just used for the Interval class. 
 * It is NOT a representation of an entire genome.
 * 
 * @author pcingola
 */
public class Genome extends Marker implements Serializable, Iterable<Chromosome> {

	private static final long serialVersionUID = -330362012383572257L;

	long length = -1;
	String species;
	String version;
	String fastaDir;
	ArrayList<String> chromosomeNames;
	ArrayList<String> chromosomeNamesSorted = null;
	String chromoFastaFiles[];
	HashMap<String, Chromosome> chromosomes;
	Genes genes; // All genes, transcripts, exons, UTRs, CDS, etc.
	Boolean codingInfo = null; // Do we have coding info from genes?

	/**
	 * Create a genome from a faidx file.
	 * See "samtools faidx" command (reference http://samtools.sourceforge.net/samtools.shtml) 
	 * 
	 * @param genomeName : Genome's name (version)
	 * @param faidxFile : FAI file used to create all chromosomes
	 * @return
	 */
	public static Genome createFromFaidx(String genomeName, String faidxFile) {
		Genome genome = new Genome(genomeName);

		// Read the whole file
		String lines[] = Gpr.readFile(faidxFile).split("\n");
		for (String line : lines) {
			String vals[] = line.split("\t");
			String chrName = vals[0];
			int len = Gpr.parseIntSafe(vals[1]);

			// Create chromo
			Chromosome chromosome = new Chromosome(genome, 0, len, 1, chrName);
			genome.add(chromosome);
		}
		return genome;

	}

	public Genome(String version) {
		super(null, Integer.MIN_VALUE, Integer.MAX_VALUE, 1, version);
		this.version = version;
		chromosomeNames = new ArrayList<String>();
		chromosomes = new HashMap<String, Chromosome>();
		genes = new Genes(this);
	}

	public Genome(String version, Properties properties) {
		super(null, Integer.MIN_VALUE, Integer.MAX_VALUE, 1, version);
		this.version = version;
		genes = new Genes(this);

		species = properties.getProperty(version + ".genome");
		if (species == null) throw new RuntimeException("Property: '" + version + ".genome' not found");
		species = species.trim();

		chromosomeNames = new ArrayList<String>();
		String[] chromosomeNames = propertyToStringArray(properties, version + ".chromosomes");

		// Fasta file & dir (optional)
		if (properties.getProperty(version + ".fasta_dir") != null) fastaDir = properties.getProperty(version + ".fasta_dir").trim();
		else fastaDir = "";
		if (properties.getProperty(version + ".chromo_fasta_files") != null) chromoFastaFiles = propertyToStringArray(properties, version + ".chromo_fasta_files");
		else chromoFastaFiles = new String[0];

		chromosomes = new HashMap<String, Chromosome>();
		for (String chName : chromosomeNames)
			add(new Chromosome(this, 0, 0, 0, chName));

	}

	public Genome(String species, String version) {
		super(null, Integer.MIN_VALUE, Integer.MAX_VALUE, 1, version);
		this.species = species;
		this.version = version;
		chromosomeNames = new ArrayList<String>();
		chromoFastaFiles = new String[0];
		chromosomes = new HashMap<String, Chromosome>();
	}

	/**
	 * Add a chromosome
	 * @param chromo
	 */
	public void add(Chromosome chromo) {
		chromosomeNames.add(chromo.getId());
		chromosomes.put(chromo.getId(), chromo);
	}

	/**
	 * Get a sorted list of chromosomes
	 * @return
	 */
	public List<String> chromosomeNamesSorted() {
		if (chromosomeNamesSorted != null) return chromosomeNamesSorted; // Already done? => return previous result

		// Sort chromosomes by name
		ArrayList<Chromosome> chromosArr = new ArrayList<Chromosome>(chromosomes.size());
		chromosArr.addAll(chromosomes.values());
		Collections.sort(chromosArr);

		// Create a list and add all names to list
		chromosomeNamesSorted = new ArrayList<String>();
		for (int i = 0; i < chromosArr.size(); i++)
			chromosomeNamesSorted.add(chromosArr.get(i).getId());

		return chromosomeNamesSorted;
	}

	public String[] getChromoFastaFiles() {
		return chromoFastaFiles;
	}

	/**
	 * Find chromosome 'chromoName'
	 * @param chromoName
	 * @return
	 */
	public Chromosome getChromosome(String chromoName) {
		String ch = Chromosome.simpleName(chromoName);
		return chromosomes.get(ch);
	}

	public String[] getChromosomeNames() {
		return chromosomeNames.toArray(new String[0]);
	}

	public String getFastaDir() {
		return fastaDir;
	}

	public Genes getGenes() {
		return genes;
	}

	/**
	 * Create a sorted list of genes (sorted by gene Id)
	 * @return
	 */
	public List<Gene> getGenesSorted() {
		ArrayList<Gene> genesSorted = new ArrayList<Gene>();
		genesSorted.addAll(genes.values());
		Collections.sort(genesSorted, new Comparator<Gene>() {

			@Override
			public int compare(Gene gene1, Gene gene2) {
				return gene1.getId().compareTo(gene2.getId());
			}
		});

		return genesSorted;
	}

	/**
	 * Create a sorted list of genes (sorted by genomic position)
	 * @return
	 */
	public List<Gene> getGenesSortedPos() {
		ArrayList<Gene> genesSorted = new ArrayList<Gene>();
		genesSorted.addAll(genes.values());
		Collections.sort(genesSorted);
		return genesSorted;
	}

	/**
	 * Get or create a chromosome
	 * @param chromoName
	 * @return
	 */
	public Chromosome getOrCreateChromosome(String chromoName) {
		String ch = Chromosome.simpleName(chromoName);
		Chromosome chr = getChromosome(ch);

		if (chr == null) {
			chr = new Chromosome(this, 0, 0, 1, chromoName);
			add(chr);
		}
		return chr;
	}

	public String getSpecies() {
		return species;
	}

	public String getVersion() {
		return version;
	}

	/**
	 * Is this chromosome in this genome?
	 * @param chromo
	 * @return
	 */
	public boolean hasChromosome(String chromo) {
		for (String ch : chromosomeNames)
			if (ch.equals(chromo)) return true;
		return false;
	}

	/**
	 * Do we have coding info from genes?
	 * @return
	 */
	public boolean hasCodingInfo() {
		// Is this already calculated?
		if (codingInfo == null) {
			int countCoding = 0;

			for (Gene gene : genes)
				if (gene.isProteinCoding()) countCoding++;

			codingInfo = (countCoding != 0);
		}

		return codingInfo;
	}

	@Override
	public Iterator<Chromosome> iterator() {
		return chromosomes.values().iterator();
	}

	/**
	 * Total genome length: add all chromosomes
	 * @return
	 */
	public long length() {
		if (length <= 0) {
			length = 0;
			for (Chromosome chr : chromosomes.values())
				length += chr.getEnd() - chr.getStart() + 1;
		}

		return length;
	}

	/**
	 * Parse a comma separated property as a string array
	 * @param properties
	 * @param attr
	 * @return
	 */
	String[] propertyToStringArray(Properties properties, String attr) {
		String value = properties.getProperty(attr);
		if (value == null) return new String[0];

		String values[] = value.split("[\\s+,]");
		LinkedList<String> list = new LinkedList<String>();
		for (String val : values)
			if (val.length() > 0) list.add(val);

		return list.toArray(new String[0]);
	}

	/**
	 * Read the whole genome sequence into memory
	 * @param fastaFile : Path to a Fasta file 
	 * @return true if it was successful
	 */
	public boolean readGenomeSequence(String fastaFile) {
		// Read fasta sequence
		FastaFileIterator ffi = new FastaFileIterator(fastaFile);
		for (String seq : ffi) {
			String chrName = ffi.getName();
			Chromosome chromosome = getChromosome(chrName);
			if (chromosome != null) {
				chromosome.setSequence(seq);
			} else {
				// Chromosome not found, create a new one
				chromosome = new Chromosome(this, 0, seq.length(), 1, chrName);
				chromosome.setSequence(seq);
				add(chromosome);
			}
		}

		return true;
	}

	/**
	 * Show number of genes, transcripts & exons 
	 * @return true : If there is an error condition (most exons do not have sequences)
	 */
	public boolean showStats() {
		int exonSeq = 0, exonNoSeq = 0;
		int countGenes = 0, countGenesProteinCoding = 0;
		int countTranscripts = 0, countTranscriptsProteinCoding = 0;
		int countExons = 0, countCds = 0;
		Genes genes = getGenes();

		for (Gene g : genes) {
			countGenes++;
			if (g.isProteinCoding()) countGenesProteinCoding++;

			for (Transcript tr : g) {
				if (tr.isProteinCoding()) countTranscriptsProteinCoding++;

				int numCds = tr.getCds().size();
				int numExons = tr.subintervals().size();

				countTranscripts++;
				countExons += numExons;
				countCds += numCds;

				for (Exon e : tr) {
					if (e.getSequence().isEmpty()) exonNoSeq++;
					else exonSeq++;
				}
			}
		}

		// Show summary
		System.out.println("# Has protein coding info    : " + hasCodingInfo());
		System.out.println("# Genes                      : " + countGenes);
		System.out.println("# Protein coding genes       : " + countGenesProteinCoding);
		System.out.println("# Transcripts                : " + countTranscripts);
		System.out.println("# Protein coding transcripts : " + countTranscriptsProteinCoding);
		System.out.println("# Cds                        : " + countCds);
		System.out.println("# Exons                      : " + countExons);
		System.out.println("# Exons with sequence        : " + exonSeq);
		System.out.println("# Exons without sequence     : " + exonNoSeq);
		return exonSeq < exonNoSeq;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder(version + ": " + species);

		if (chromosomeNames.size() > 0) {
			sb.append("\n\tChromosomes: ");
			for (String chr : chromosomeNames)
				sb.append(chr + " ");
			sb.append("\n");
		}

		if ((chromoFastaFiles != null) && (chromoFastaFiles.length > 0)) {
			sb.append("\tFasta files: ");
			for (String ff : chromoFastaFiles)
				sb.append(fastaDir + "/" + ff + " ");

			sb.append("\n");
		}

		return sb.toString();
	}

}