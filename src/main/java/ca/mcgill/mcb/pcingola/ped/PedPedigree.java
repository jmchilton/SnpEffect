package ca.mcgill.mcb.pcingola.ped;

import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;

import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * A pedigree of PedEntries
 * 
 * @author pcingola
 */
public class PedPedigree implements Iterable<TfamEntry> {

	boolean verbose = false;
	HashMap<String, TfamEntry> tfamById = new HashMap<String, TfamEntry>();
	PlinkMap plinkMap;

	public PedPedigree() {
		tfamById = new HashMap<String, TfamEntry>();
	}

	public PedPedigree(String tfamFileName) {
		tfamById = new HashMap<String, TfamEntry>();
		loadTfam(tfamFileName);
	}

	/**
	 * Add an entry to this pedigree
	 * @param tfamEntry
	 */
	public void add(TfamEntry tfamEntry) {
		tfamById.put(tfamEntry.getId(), tfamEntry);
	}

	public TfamEntry get(String id) {
		return tfamById.get(id);
	}

	public PlinkMap getPlinkMap() {
		return plinkMap;
	}

	@Override
	public Iterator<TfamEntry> iterator() {
		return tfamById.values().iterator();
	}

	public Set<String> keySet() {
		return tfamById.keySet();
	}

	/**
	 * Load a pedigree from a PED and MAP file pair
	 * 
	 * @param pedFileName
	 */
	public void load(String pedFileName) {
		String pedBaseFileName = Gpr.removeExt(pedFileName);
		String mapFile = pedBaseFileName + ".map";

		PedFileIterator pedFile = new PedFileIterator(pedFileName, mapFile);

		// Load all entries for this family
		int count = 1;
		for (PedEntry pe : pedFile) {
			if (verbose) Gpr.showMarkStderr(count++, 1);
			add(pe);
		}

		plinkMap = pedFile.getPlinkMap();
	}

	/**
	 * Load a TFAM file
	 * @param tfamFileName
	 */
	public void loadTfam(String tfamFileName) {
		String tfamStr = Gpr.readFile(tfamFileName);
		if (tfamStr.isEmpty()) throw new RuntimeException("Cannot read file '" + tfamFileName + "'");

		for (String line : tfamStr.split("\n")) {
			TfamEntry te = new TfamEntry(line);
			add(te);
		}
	}

	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}

	public int size() {
		return tfamById.size();
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (TfamEntry te : this)
			sb.append(te.toString() + "\n");
		return sb.toString();
	}

	public Collection<TfamEntry> values() {
		return tfamById.values();
	}
}
