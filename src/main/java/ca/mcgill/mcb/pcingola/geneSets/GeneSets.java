package ca.mcgill.mcb.pcingola.geneSets;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import ca.mcgill.mcb.pcingola.geneOntology.GoTerm;
import ca.mcgill.mcb.pcingola.geneOntology.GoTerms;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * A collection of GeneSets
 * 
 * @author Pablo Cingolani
 */
@SuppressWarnings("serial")
public class GeneSets implements Iterable<GeneSet>, Serializable {

	public static boolean debug = false; // Debug mode for this class?
	public static double LOG2 = Math.log(2); // We use this constant often
	public static long PRINT_SOMETHING_TIME = 5000; // Print something every X seconds
	static int warnCount = 0;

	boolean verbose = false; // Verbose mode for this class?
	boolean doNotAddIfNotInGeneSet = false; // Do not add genes that don't belong to geneset
	int maxRank; // Maximum rank in this collection
	String label; // Label for this set of nodes 
	HashSet<String> genes; // All genes in this experiment
	HashMap<String, GeneSet> geneSetsByName; // Gene sets indexed by GeneSet.name
	HashMap<String, HashSet<GeneSet>> geneSetsByGene; // Gene sets indexed by gene name
	HashSet<String> interestingGenes; // Interesting genes in this experiment
	HashMap<String, Integer> rankByGene; // Ranked genes
	HashMap<String, Double> valueByGene;

	/**
	 * Create gene sets form GoTerms
	 * @param goTerms : GoTerms to use
	 */
	public static GeneSets factory(GoTerms goTerms) {
		GeneSets geneSets = new GeneSets();

		for (GoTerm gt : goTerms) {
			// Create gene set
			GeneSet geneSet = new GeneSet(gt.getAcc(), gt.getDescription(), geneSets);

			// Add all genes
			for (String id : gt)
				geneSet.addGene(id);

			// Add to Gene Sets
			geneSets.add(geneSet);
		}

		return geneSets;
	}

	/**
	 * Default constructor
	 */
	public GeneSets() {
		geneSetsByName = new HashMap<String, GeneSet>();
		interestingGenes = new HashSet<String>();
		genes = new HashSet<String>();
		rankByGene = new HashMap<String, Integer>();
		geneSetsByGene = new HashMap<String, HashSet<GeneSet>>();
		maxRank = 0;
	}

	public GeneSets(String msigDb) {
		geneSetsByName = new HashMap<String, GeneSet>();
		interestingGenes = new HashSet<String>();
		genes = new HashSet<String>();
		rankByGene = new HashMap<String, Integer>();
		geneSetsByGene = new HashMap<String, HashSet<GeneSet>>();
		maxRank = 0;
		loadMSigDb(msigDb, false);
	}

	/**
	 * Add a gene set
	 * @param geneSetName
	 * @param geneSet
	 */
	public void add(GeneSet geneSet) {
		// Add
		geneSetsByName.put(geneSet.getName().toUpperCase(), geneSet);

		// Add all genes
		for (String gene : geneSet)
			add(gene, geneSet);

		geneSet.setGeneSets(this);
	}

	/**
	 * Add a gene and aliases
	 */
	public boolean add(String gene) {
		return genes.add(gene);
	}

	/**
	 * Add a gene and it's corresponding gene set
	 * @param gene
	 * @param geneSet
	 * @return
	 */
	public boolean add(String gene, GeneSet geneSet) {
		HashSet<GeneSet> listgs = geneSetsByGene.get(gene);
		if (listgs == null) {
			listgs = new HashSet<GeneSet>();
			geneSetsByGene.put(gene, listgs);
		}
		listgs.add(geneSet);

		return genes.add(gene);
	}

	/**
	 * Add a symbol as 'interesting' gene (to every corresponding GeneSet in this collection) 
	 * @param gene : symbol's ID
	 * @param rank : symbol's rank
	 * @returns : true if it was added OK, false on error.
	 */
	public boolean add(String gene, int rank) {
		boolean ok = true;

		// Sanity check
		if (!genes.contains(gene)) {
			if (verbose) System.err.println("WARNING: Trying to add ranked gene. Gene  '" + gene + "' does not exist in GeneSets! " + (doNotAddIfNotInGeneSet ? "Ignored." : "Added anyway."));
			ok = false;
			if (doNotAddIfNotInGeneSet) return ok;
		}

		rankByGene.put(gene, rank); // Add gene -> rank pair 
		interestingGenes.add(gene);
		if (maxRank < rank) maxRank = rank;

		return ok;
	}

	/**
	 * Checks that every symboolID is in the set (as 'interesting' genes)
	 * @param intGenes : A set of interesting genes
	 * Throws an exception on error
	 */
	public void checkInterestingGenes(Set<String> intGenes) {
		if (debug) Timer.showStdErr("Checking genes (" + intGenes.size() + ") : " + intGenes);

		// Check that genes contains interestingGenes
		if (!intGenes.containsAll(interestingGenes)) { throw new RuntimeException("Not every gene in :" + label + " as an interesting symbol"); }

		// Check that every interesting symbol in DAG is from genes
		if (!interestingGenes.containsAll(intGenes)) { throw new RuntimeException("Not every gene marked as interesting in " + label + " is from intGenes\n\tInteresting genes(" + interestingGenes.size() + "): " + interestingGenes + "\n\tintGenes(" + intGenes.size() + "): " + intGenes); }

		// Are ranks being used? => Check them
		if ((rankByGene != null) && (rankByGene.keySet().size() > 0)) {
			int maxRankTmp = intGenes.size();

			// Check that every rank is being used
			int ranksUsed[] = new int[maxRankTmp + 1];
			for (int i = 0; i < maxRankTmp; i++)
				ranksUsed[i] = 0;

			// Check that every interestingSymbolId is ranked (if ranks are being used)
			for (String gene : intGenes) {
				Integer rank = rankByGene.get(gene);
				if ((rank == null) || (rank <= 0) || (rank > maxRankTmp)) { throw new RuntimeException("Invalid rank for gene:" + gene + ", rank:" + rank + "(should be [1," + maxRankTmp + "]"); }
				ranksUsed[rank]++;
			}

			for (int rank = 1; rank < maxRankTmp; rank++) {
				if (ranksUsed[rank] != 1) { throw new RuntimeException("Rank number " + rank + " is used " + ranksUsed[rank] + " times (should be used exactly 1 time)"); }
			}
		}
	}

	/**
	 * Produce a GeneSet based on a list of GeneSets and a 'mask'
	 * 
	 * @param geneSetList : A list of GeneSets
	 * @param activeSets : An integer (binary mask) that specifies weather a set in the list should be taken into account or not. The operation performed is: 
	 * 
	 * 		Intersection{ GeneSets where mask_bit == 1 } - Union{ GeneSets where mask_bit == 0 } )
	 * 
	 * where the minus sign '-' is actually a 'set minus' operation. This operation is done for both sets
	 * in GeneSet (i.e. genes and interestingGenes)
	 * 
	 * @return A GeneSet
	 */
	public GeneSet disjointSet(List<GeneSet> geneSetList, int activeSets) {
		//---
		// Produce intersections (for each term in the list)
		//---
		GeneSet gtUnion = new GeneSet("UNION", "UNION", null);
		GeneSet gtIntersect = new GeneSet("INTERSECTION", "INTERSECTION", null);

		int i = 0;
		boolean firstIntersection = true;
		for (GeneSet geneSet : geneSetList) {
			// Extract the i_th bit from 'activeSets'
			boolean biti = (activeSets & (1L << i)) > 0;

			if (biti) { // Is this bit is 1? => perform an intersection

				if (firstIntersection) { // Initialize intersection set (otherwise all intersections are empty)
					gtIntersect.union(geneSet);
					firstIntersection = false;
				} else {
					gtIntersect.intersection(geneSet);
					// Are we done? (if the intersection set is empty, it doesn't make any sense to continue
					if (gtIntersect.getGeneCount() <= 0) return gtIntersect;
				}
			} else gtUnion.union(geneSet);

			i++;
		}

		// Now extract the 'union' set from the intersection set (i.e. perform a 'set minus' operation)
		gtIntersect.setMinus(gtUnion);

		return gtIntersect;
	}

	/**
	 * Iterate through each GeneSet in this GeneSets
	 */
	public List<GeneSet> geneSetsSorted() {
		LinkedList<GeneSet> ll = new LinkedList<GeneSet>(geneSetsByName.values());
		Collections.sort(ll);
		return ll;
	}

	/**
	 * Gene sets sorted by size (if same size, sort by name).
	 * @param reverse : Reverse size sorting (does not affect name sorting)
	 * @return
	 */
	public List<GeneSet> geneSetsSortedSize(final boolean reverse) {
		ArrayList<GeneSet> ll = new ArrayList<GeneSet>(geneSetsByName.values());
		Collections.sort(ll, new Comparator<GeneSet>() {

			@Override
			public int compare(GeneSet gs1, GeneSet gs2) {
				// Compare by size
				int diff = gs1.size() - gs2.size();
				if (diff != 0) return (reverse ? -diff : diff);

				// Same size? Sozr by name
				return gs1.getName().compareTo(gs2.getName());
			}
		});
		return ll;
	}

	/**
	 * How many genes do we have?
	 * @return
	 */
	public int getGeneCount() {
		if (genes == null) { return 0; }
		return genes.size();
	}

	/** 
	 * Get all genes in this set
	 * @return
	 */
	public Set<String> getGenes() {
		return genes;
	}

	/**
	 * Get a gene set named 'geneSetName'
	 * @param geneSetName
	 * @return
	 */
	public GeneSet getGeneSet(String geneSetName) {
		return geneSetsByName.get(geneSetName.toUpperCase());
	}

	public int getGeneSetCount() {
		if (geneSetsByName == null) { return 0; }
		return geneSetsByName.size();
	}

	/**
	 * All gene sets that this gene belongs to
	 * @param gene
	 * @return
	 */
	public HashSet<GeneSet> getGeneSetsByGene(String gene) {
		return geneSetsByGene.get(gene);
	}

	public HashMap<String, GeneSet> getGeneSetsByName() {
		return geneSetsByName;
	}

	public HashSet<String> getInterestingGenes() {
		return interestingGenes;
	}

	public int getInterestingGenesCount() {
		return interestingGenes.size();
	}

	public String getLabel() {
		return label;
	}

	public int getMaxRank() {
		// Calculate it if needed
		if (maxRank <= 0) { // Find max rank used 

			for (String gene : rankByGene.keySet()) {
				int rank = rankByGene.get(gene);
				if (rank > maxRank) maxRank = rank;
			}
		}
		return maxRank;
	}

	/**
	 * Get gene's rank
	 * @param gene
	 * @return
	 */
	public int getRank(String gene) {
		Integer rank = rankByGene.get(gene);
		if (rank == null) { return 0; }
		return rank;
	}

	public HashMap<String, Integer> getRankByGene() {
		return rankByGene;
	}

	/**
	 * How many gene sets have ranked genes (i.e. rank sum > 0)
	 * @return Number of gene set such that rankSum > 0
	 */
	public int getRankedSetsCount() {
		int count = 0;
		for (GeneSet gs : this)
			if (gs.rankSum() > 0) count++;
		return count;
	}

	/**
	 * Get experimental value
	 * @param gene
	 * @return
	 */
	public double getValue(String gene) {
		Double val = valueByGene.get(gene);
		if (val == null) { return 0; }
		return val;
	}

	public HashMap<String, Double> getValueByGene() {
		return valueByGene;
	}

	/**
	 * Iterate through each GeneSet in this GeneSets
	 */
	@Override
	public Iterator<GeneSet> iterator() {
		return geneSetsByName.values().iterator();
	}

	/**
	 * Iterate through each GeneSet in this GeneSets
	 */
	public Iterator<GeneSet> iteratorSorted() {
		return geneSetsSorted().iterator();

	}

	public Set<String> keySet() {
		return geneSetsByName.keySet();
	}

	/**
	 * Select a number of GeneSets
	 * @param numberToSelect
	 * @return
	 */
	public List<GeneSet> listTopTerms(int numberToSelect) {
		LinkedList<GeneSet> list = new LinkedList<GeneSet>();

		// Create a list of terms (to be ordered)
		int i = 0;
		LinkedList<GeneSet> ll = new LinkedList<GeneSet>();
		for (String geneSetName : keySet())
			ll.add(getGeneSet(geneSetName));

		Collections.sort(ll);

		for (GeneSet geneSet : ll)
			if (i++ < numberToSelect) list.add(geneSet);

		return list;
	}

	/**
	 * Reads a file with a list of genes and experimental values.
	 * Format: "gene \t value \n"
	 * @param fileName 
	 * @return A list of genes not found
	 */
	@SuppressWarnings("unchecked")
	public List<String> loadExperimentalValues(String fileName, boolean maskException) {
		LinkedList<String> notFound = new LinkedList<String>();

		if (verbose) System.err.println("Reading 'ranked' genes from file: '" + fileName + "'");

		// First: Initialize
		reset();

		//---
		// Read genes and values
		//---

		// Read file
		try {
			// Open file and initialize buffers
			BufferedReader inFile = new BufferedReader(new FileReader(fileName));
			String line;

			// Read each line and add it to 'valueByGene'
			while ((line = inFile.readLine()) != null) {
				// line = line.trim();
				if (!line.startsWith("#")) { // Skip comments

					String fields[] = line.split("\t");
					String gene = fields[0];
					if (!gene.isEmpty()) {
						Double value = Double.parseDouble(fields[1]);
						valueByGene.put(gene, value);

						// Is this gene in this gene set collection
						if (!genes.contains(gene)) notFound.add(gene);
					}
				}
			}

			// OK, finished
			inFile.close();
		} catch (IOException e) {
			if (maskException) { return null; }
			throw new RuntimeException(e);
		}

		//---
		// Sort by experimental value then rank them
		//---
		LinkedList<String> geneNames = new LinkedList<String>(valueByGene.keySet());
		Collections.sort(geneNames, new CompareByValue(valueByGene, false));
		int rank = 1, errorsRank = 0;
		for (String gene : geneNames)
			if (!add(gene, rank++)) errorsRank++;

		if (verbose) System.err.println("Genes added: " + geneNames.size() + "\nErrors adding ranked / interesting genes: " + errorsRank);
		return notFound;
	}

	/**
	 * Read an MSigDBfile and add every Gene set (do not add relationships between nodes in DAG)
	 * @param gmtFile
	 * @param geneSetType
	 */
	public boolean loadMSigDb(String gmtFile, boolean maskException) {
		try {
			if (verbose) System.err.println("Reading gene sets file: '" + gmtFile + "'");

			genes = new HashSet<String>(); // Reset genes
			geneSetsByName = new HashMap<String, GeneSet>(); // and genesets by name

			// Open file and initialize buffers
			BufferedReader inFile = new BufferedReader(new FileReader(gmtFile));
			String line;
			int lineNum;

			// Read each line
			for (lineNum = 0; (line = inFile.readLine()) != null; lineNum++) {
				line = line.trim();
				if (!line.startsWith("#")) { // Skip comments

					String fields[] = line.split("\t");

					// Parse name & description
					String geneSetName = fields[0];
					String description = fields[1];

					// Sanity check: Does gene set already exist?
					if (getGeneSet(geneSetName) != null) Gpr.debug("Error: File '" + gmtFile + "' line " + lineNum + ". Gene set name '" + geneSetName + "' duplicated.");

					// Create geneSet and all genes
					GeneSet gs = new GeneSet(geneSetName, description, this);
					for (int i = 2; i < fields.length; i++)
						gs.addGene(fields[i]);

					// Add gene set to collection
					add(gs);
				}
			}

			// OK, finished
			inFile.close();

			if (verbose) System.err.println("GeneSets added: " + lineNum);

		} catch (IOException e) {
			if (maskException) return false;
			throw new RuntimeException(e);
		}

		return true;
	}

	/**
	 * Remove a GeneSet
	 */
	public void removeGeneSet(String geneSetName) {
		geneSetsByName.remove(geneSetName);
	}

	/**
	 * Reset every 'interesting' gene or ranked gene (on every single GeneSet in this GeneSets)
	 */
	public void reset() {
		interestingGenes = new HashSet<String>();
		rankByGene = new HashMap<String, Integer>();
		valueByGene = new HashMap<String, Double>();
		maxRank = 0;

		for (GeneSet gt : this) {
			gt.reset();
		}
	}

	/**
	 * Save gene sets file for GSEA analysis
	 * Format specification: http://www.broad.mit.edu/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29
	 * 
	 * @param fileName
	 */
	public void saveGseaGeneSets(String fileName) {
		// Create a string with all the data
		StringBuffer out = new StringBuffer();
		for (GeneSet gt : this) // Save GeneSet that have at least 1 gene
		{
			if (gt.getGenes().size() > 0) {
				out.append(gt.getName() + "\t" + gt.getName() + "\t");

				// Add all genes for this GeneSet
				for (String gene : gt.getGenes())
					out.append("gene_" + gene + "\t");
				out.append("\n");
			}

			// Save it 
		}
		Gpr.toFile(fileName, out);
	}

	public void setDoNotAddIfNotInGeneSet(boolean doNotAddIfNotInGeneSet) {
		this.doNotAddIfNotInGeneSet = doNotAddIfNotInGeneSet;
	}

	public void setGeneSetByName(HashMap<String, GeneSet> geneSets) {
		geneSetsByName = geneSets;
	}

	public void setInterestingGenes(HashSet<String> interestingGenesIdSet) {
		interestingGenes = interestingGenesIdSet;
	}

	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}

	/**
	 * Sort and rank based on experimental value
	 * @param orderAscending
	 */
	@SuppressWarnings("unchecked")
	public void sortAndRank(HashMap<String, Double> valueByGene, boolean orderAscending) {
		if (verbose) System.err.println("Sorting: " + (orderAscending ? "ascentding" : "descending"));

		// Sort genes, rank them and add them to GeneSets
		LinkedList<String> geneNames = new LinkedList<String>(valueByGene.keySet());
		Collections.sort(geneNames, new CompareByValue(valueByGene, orderAscending));
		int rank = 1;
		int showEvery = geneNames.size() / 100 + 1;
		for (String gene : geneNames) {
			if (verbose && (((rank - 1) % showEvery) == 0)) System.err.println("\tRank: " + rank + "\tGene: " + gene + "\tValue: " + valueByGene.get(gene));
			if (add(gene, rank)) rank++;
		}

		if (verbose) System.err.println("Total genes added: :" + interestingGenes.size() + "\tMax rank: " + getMaxRank());
	}

	@Override
	public String toString() {
		StringBuffer sb = new StringBuffer();

		for (GeneSet gs : geneSetsSorted())
			sb.append(gs.toStringAll() + "\n");

		return sb.toString();
	}

	public Collection<GeneSet> values() {
		return geneSetsByName.values();
	}
}
