package ca.mcgill.mcb.pcingola.interval;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.zip.GZIPOutputStream;

import ca.mcgill.mcb.pcingola.fileIterator.LineFileIterator;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.SnpEffectPredictor;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * Serialize markers to (and from) file
 * 
 * @author pcingola
 */
public class MarkerSerializer {

	PrintStream outFile;
	int lineNum;
	String line;
	int parsedField;
	String fields[];
	int currId = 0;

	HashMap<Integer, Marker> markerById;
	HashMap<Marker, Integer> idByMarker;

	public static void main(String[] args) throws IOException {
		String fileName = Gpr.HOME + "/zzz.txt.gz";

		// Read database
		Timer.showStdErr("Loading");
		Config config = new Config("testHg3763ChrY");
		SnpEffectPredictor sep = SnpEffectPredictor.load(config);

		Markers markers = new Markers();
		markers.add(sep.getGenome());

		for (Chromosome chr : sep.getGenome())
			markers.add(chr);

		for (Gene g : sep.getGenome().getGenes())
			markers.add(g);

		// Write
		MarkerSerializer is = new MarkerSerializer();
		Timer.showStdErr("Writing to " + fileName);
		is.save(fileName, markers);

		// Read
		Timer.showStdErr("Reading to " + fileName);
		is.load(fileName);
		Timer.showStdErr("Done.");

	}

	public MarkerSerializer() {
		markerById = new HashMap<Integer, Marker>();
		idByMarker = new HashMap<Marker, Integer>();
	}

	public int getIdByMarker(Marker m) {
		return idByMarker.get(m);
	}

	public Marker getMarkerById(int id) {
		return markerById.get(id);
	}

	public String getNextField() {
		return fields[parsedField++];
	}

	public boolean getNextFieldBoolean() {
		return Gpr.parseBoolSafe(getNextField());
	}

	public int getNextFieldInt() {
		return Gpr.parseIntSafe(getNextField());
	}

	public Marker getNextFieldMarker() {
		return getMarkerById(getNextFieldInt());
	}

	public Markers getNextFieldMarkers() {
		Markers markers = new Markers();
		String fieldIdsStr = getNextField();
		if (fieldIdsStr.isEmpty()) return markers;

		String fieldIds[] = fieldIdsStr.split(",");
		for (String idStr : fieldIds) {
			int id = Gpr.parseIntSafe(idStr);
			Marker m = getMarkerById(id);
			if (m != null) markers.add(m);
			else throw new RuntimeException("Marker '" + id + "' not found. This should never happen!");
		}
		return markers;
	}

	public int getNextId() {
		return ++currId;
	}

	/**
	 * Load data from file
	 * @param fileName
	 */
	public void load(String fileName) {
		//---
		// Load data from file
		//---
		LineFileIterator lfi = new LineFileIterator(fileName);
		int lineNum = 0;
		for (String l : lfi) {
			line = l;
			parsedField = 0;
			fields = line.split("\t", -1);

			// Parse field type
			String typeStr = fields[0];
			EffectType type = EffectType.valueOf(typeStr);

			// Parse serialization id
			String idStr = fields[1];
			int id = Gpr.parseIntSafe(idStr);

			Marker m = null;
			switch (type) {
			case GENE:
				m = new Gene();
				break;
			case UTR_3_PRIME:
				m = new Utr3prime();
				break;
			case UTR_5_PRIME:
				m = new Utr5prime();
				break;
			case EXON:
				m = new Exon();
				break;
			case CDS:
				m = new Cds();
				break;
			case TRANSCRIPT:
				m = new Transcript();
				break;
			case CHROMOSOME:
				m = new Chromosome();
				break;
			case GENOME:
				m = new Genome();
				break;
			default:
				throw new RuntimeException("Unimplemented for type '" + type + "'");
			}

			try {
				// Parse line
				m.serializeParse(this);
				System.out.println(m.getId());
			} catch (Throwable t) {
				t.printStackTrace();
				throw new RuntimeException("Error parsing line " + (lineNum + 1) + " from file '" + fileName + "'\n\t" + line + "\n\tField [" + parsedField + "] : '" + (parsedField < fields.length ? fields[parsedField] : "-") + "'", t);
			}

			// Add to hash
			markerById.put(id, m);
			lineNum++;
		}

		//--- 
		// Assign parents
		//---
		for (Marker m : markerById.values()) {
			MarkerParentId mpid = (MarkerParentId) m.getParent();
			int parentId = mpid.getParentId();
			Marker parent = getMarkerById(parentId);
			m.setParent(parent);
		}
	}

	public String save(Iterable<Marker> markersCol) {
		StringBuilder idStr = new StringBuilder();
		for (Marker m : markersCol) {
			int id = save(m);
			if (idStr.length() > 0) idStr.append(",");
			idStr.append(id);
		}
		return idStr.toString();
	}

	/**
	 * Save a marker
	 * @param m
	 */
	public int save(Marker m) {
		if (m == null) return -1;
		if (idByMarker.containsKey(m)) return idByMarker.get(m); // Already done

		int id = getNextId();
		idByMarker.put(m, id);
		String line = m.serializeSave(this);
		outFile.print(line + "\n");
		lineNum++;
		return id;
	}

	/**
	 * Save a genome
	 * @param fileName
	 * @param genome
	 */
	public void save(String fileName, Genome genome) {
		Markers markers = new Markers();
		markers.add(genome);

		for (Chromosome chr : genome)
			markers.add(chr);

		for (Gene g : genome.getGenes())
			markers.add(g);

		// Write
		save(fileName, markers);
	}

	/**
	 * Save data to file
	 * @param fileName
	 */
	public void save(String fileName, Markers markers) {
		try {
			lineNum = 0;
			currId = 0;
			outFile = new PrintStream(new GZIPOutputStream(new FileOutputStream(fileName)));
			for (Marker m : markers)
				save(m);
			outFile.close();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
}
