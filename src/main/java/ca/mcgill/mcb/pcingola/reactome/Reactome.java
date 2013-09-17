package ca.mcgill.mcb.pcingola.reactome;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import ca.mcgill.mcb.pcingola.collections.AutoHashMap;
import ca.mcgill.mcb.pcingola.fileIterator.LineFileIterator;
import ca.mcgill.mcb.pcingola.gtex.Gtex;
import ca.mcgill.mcb.pcingola.gtex.GtexExperiment;
import ca.mcgill.mcb.pcingola.reactome.events.BlackBoxEvent;
import ca.mcgill.mcb.pcingola.reactome.events.CatalystActivity;
import ca.mcgill.mcb.pcingola.reactome.events.Complex;
import ca.mcgill.mcb.pcingola.reactome.events.Depolymerisation;
import ca.mcgill.mcb.pcingola.reactome.events.Event;
import ca.mcgill.mcb.pcingola.reactome.events.Pathway;
import ca.mcgill.mcb.pcingola.reactome.events.Polymerisation;
import ca.mcgill.mcb.pcingola.reactome.events.Reaction;
import ca.mcgill.mcb.pcingola.stats.CountByType;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * Load reactome data from TXT files
 * 
 * @author pcingola
 *
 */
public class Reactome implements Iterable<Entity> {

	public static final int SHOW_EVERY = 10000;
	public static final double EPSILON = 1E-6;
	public static final double MAX_CONVERGENCE_DIFFERENCE = 1E-3;
	public static final int MAX_ITERATIONS = 1000;

	String dirName;

	HashMap<String, Entity> entityById;

	HashMap<String, String> objectType;

	HashMap<String, String> objectName;
	AutoHashMap<String, ArrayList<Entity>> entitiesByGeneId;
	HashSet<String> entitiesGeneId = new HashSet<String>();
	Monitor monitor; // Monitor all nodes in the circuit
	Monitor monitorTrace; // Monitor a specific set of nodes (usually one node and all it's predecesors)

	/**
	 * Find and return first Regexp occurrence (null if nothing is found)
	 * @param pattern
	 * @param str
	 * @return
	 */
	static String findRegexp(Pattern pattern, String str) {
		Matcher m = pattern.matcher(str);
		if (!m.find()) return null;
		return m.group();
	}

	static boolean hasRegexp(Pattern pattern, String str) {
		Matcher m = pattern.matcher(str);
		return m.find();
	}

	/**
	 * Main
	 * @param args
	 */
	public static void main(String[] args) {

		//---
		// Load reactome data
		//---
		String reactomeDir = Gpr.HOME + "/snpEff/db/reactome/txt/";
		String geneIdsFile = Gpr.HOME + "/snpEff/db/reactome/gene_ids/biomart_query_uniq.txt";
		Reactome reactome = new Reactome(reactomeDir);
		reactome.load();
		reactome.loadGeneIds(geneIdsFile); // Load Gene IDs data

		//---
		// Load GTEX data
		//---
		Timer.showStdErr("Loading GTEx data");
		String gtexDir = Gpr.HOME + "/snpEff/db/GTEx";
		String gtexSamples = gtexDir + "/GTEx_Analysis_Annotations_Sample_DS__Pilot_2013_01_31.txt";
		String gtexData = gtexDir + "/gtex_norm.zzz.txt";
		Gtex gtex = new Gtex(gtexSamples, gtexData);
		gtex.getClass();

		//---
		// Simulate...!?
		//---
		for (GtexExperiment gtexExperiment : gtex) {
			if (gtexExperiment.size() > 0) {
				reactome.zzz(gtexExperiment);
			}
		}
	}

	public Reactome(String dirName) {
		this.dirName = dirName;
		entityById = new HashMap<String, Entity>();
		entitiesByGeneId = new AutoHashMap<String, ArrayList<Entity>>(new ArrayList<Entity>());
	}

	/**
	 * Add an entity <-> geneId
	 * @param entity
	 * @param geneId
	 */
	public void add(Entity entity, String geneId) {
		String key = geneId + "\t" + entity.getId();
		if (entitiesGeneId.contains(key)) return; // Don't add more then once

		entity.addGeneId(geneId);
		entitiesByGeneId.getOrCreate(geneId).add(entity);
		entitiesGeneId.add(key);
	}

	/**
	 * Iterate network until convergence
	 * 
	 * @param gtexExperiment
	 */
	boolean calc(GtexExperiment gtexExperiment) {
		boolean changed = true;
		int iteration;
		System.out.print(gtexExperiment.getTissueTypeDetail() + "\t");
		for (iteration = 0; changed && iteration < MAX_ITERATIONS; iteration++) {
			changed = false;
			HashSet<Entity> done = new HashSet<Entity>();

			for (Entity e : this) {
				double outPrev = e.getOutput();
				double out = e.calc(done);

				// Output changed?
				if (Math.abs(outPrev - out) > MAX_CONVERGENCE_DIFFERENCE) changed = true;
			}
			System.out.print(".");
		}
		System.out.println(" " + iteration);

		return changed;
	}

	/** 
	 * Create a monitor for all nodes in the circuit
	 */
	Monitor createMonitor() {
		Monitor monitor = new Monitor();
		for (Entity e : this) {
			if (!e.isFixed() && e.isReaction()) monitor.add(e);
		}

		monitor.sort();
		Gpr.debug("Monitor size: " + monitor.size());

		return monitor;
	}

	/** 
	 * Create a monitor for a subset of nodes that explain "target's" output
	 */
	Monitor createMonitor(String targetNodeId) {
		// Perform 1 iteration to get a set of all nodes required for target's output
		reset();
		Entity target = entityById.get(targetNodeId);
		HashSet<Entity> done = new HashSet<Entity>();
		target.calc(done);

		Monitor monitor = new Monitor();
		for (Entity e : done)
			monitor.add(e);
		monitor.sort();
		Gpr.debug("Monitor size: " + monitor.size());

		return monitor;
	}

	/**
	 * Get or create a new entity
	 * @param id
	 * @return
	 */
	Entity getEntity(String id) {
		// Get from hash
		Entity e = entityById.get(id);
		if (e != null) return e;

		// Not available? Create entity
		String type = objectType.get(id);
		if (type == null) throw new RuntimeException("Cannot find entity type for ID '" + id + "'");
		String name = objectName.get(id);
		int idNum = Gpr.parseIntSafe(id);

		// Create according to object type
		if (type.equals("Complex")) e = new Complex(idNum, name);
		else if (type.equals("EntityCompartment") || type.equals("Compartment") || type.equals("GO_CellularComponent")) e = new Compartment(idNum, name);
		else if (type.equals("Reaction")) e = new Reaction(idNum, name);
		else if (type.equals("BlackBoxEvent")) e = new BlackBoxEvent(idNum, name);
		else if (type.equals("Pathway")) e = new Pathway(idNum, name);
		else if (type.equals("Depolymerisation")) e = new Depolymerisation(idNum, name);
		else if (type.equals("Polymerisation")) e = new Polymerisation(idNum, name);
		else if (type.equals("CatalystActivity")) e = new CatalystActivity(idNum, name);
		else e = new Entity(idNum, name);

		// Add to maps
		entityById.put(id, e);

		return e;
	}

	@Override
	public Iterator<Entity> iterator() {
		return entityById.values().iterator();
	}

	public void load() {
		Timer.showStdErr("Loading Reactome data from directory '" + dirName + "'");

		loadDatabaseObjects(); // Load a map of all object names and types
		loadComplex2HasComponent(); // Load Complex_2_hasComponent
		loadPhysicalEntity2Compartment(); // Load compartment information
		loadPathway2HasEvent(); // Load pathway data

		loadReactionlikeEvent2Input(); // Load reaction inputs
		loadReactionlikeEvent2Output(); // Load reaction outputs
		loadReactionlikeEvent2CatalystActivity(); // Load reaction catalysts
		loadRegulation(); // Load reaction regulators

		loadCatalystActivity(); // Load catalyst

		// Remove cached data, we don't need it any more
		objectType = null;
		objectName = null;

		Timer.showStdErr("Loading finished");
	}

	/**
	 * Load catalyst activity to molecule mapping
	 */
	protected void loadCatalystActivity() {
		String name = "CatalystActivity";
		String fileName = dirName + name + ".txt";
		Timer.showStdErr("Loading " + name + " from '" + fileName + "'");

		int i = 1;
		for (String line : Gpr.readFile(fileName).split("\n")) {
			// Parse line
			String rec[] = line.split("\t");
			String id = rec[0];
			String entityId = rec[1];

			if (id.equals("DB_ID")) continue; // Skip title

			// Add event to pathway
			CatalystActivity reaction = (CatalystActivity) entityById.get(id);
			if (reaction == null) continue; // Reaction not found? Skip

			Entity e = getEntity(entityId);
			reaction.addInput(e);

			Gpr.showMark(i++, SHOW_EVERY);
		}

		System.err.println("");
		Timer.showStdErr("Total catalyst entities assigned: " + (i - 1));
	}

	/**
	 * Load complexes
	 * @param name
	 * @param fileName
	 * @param map
	 */
	protected void loadComplex2HasComponent() {
		String fileName = dirName + "Complex_2_hasComponent.txt";
		Timer.showStdErr("Loading Complex_2_hasComponent from '" + fileName + "'");

		int i = 1;
		for (String line : Gpr.readFile(fileName).split("\n")) {
			// Parse line
			String rec[] = line.split("\t");
			String id = rec[0];
			int idNum = Gpr.parseIntSafe(id);
			String componentId = rec[1];

			if (idNum == 0) continue; // Skip title

			// Get complex and add entity
			Complex c = (Complex) getEntity(id);
			c.add(getEntity(componentId));

			Gpr.showMark(i++, SHOW_EVERY);
		}

		System.err.println("");
		Timer.showStdErr("Total entities added: " + entityById.size());
	}

	/**
	 * Load objects table (populate objectType and objectName maps)
	 */
	protected void loadDatabaseObjects() {
		String fileName = dirName + "DatabaseObject.txt";
		Timer.showStdErr("Loading objects from '" + fileName + "'");

		// Ensure capacity
		int numObjects = Gpr.countLines(fileName);
		Timer.showStdErr("Counting lines from '" + fileName + "'. Total lines: " + numObjects);
		objectType = new HashMap<String, String>(numObjects);
		objectName = new HashMap<String, String>(numObjects);

		// Load objects
		int i = 1;
		LineFileIterator lfi = new LineFileIterator(fileName);
		for (String line : lfi) {
			// Parse line
			String recs[] = line.split("\t");
			if (recs.length < 3) continue;
			String id = recs[0];
			String objType = recs[1];
			String objName = recs[2];

			// Add info
			objectType.put(id, objType);
			objectName.put(id, objName);

			Gpr.showMark(i++, SHOW_EVERY);
		}

		System.err.println("");
		Timer.showStdErr("Total objects loaded: " + objectName.size());
	}

	/**
	 * Load Gene IDs data, then map geneIDs <-> Entities
	 * 
	 * @param geneIdsFile
	 */
	public void loadGeneIds(String geneIdsFile) {
		//---
		// Load Gene IDs data
		//---
		Timer.showStdErr("Loading Gene IDs from " + geneIdsFile);
		GeneIds geneIds = new GeneIds(geneIdsFile);

		// Assign to geneIDs
		Pattern patternEnsg = Pattern.compile("ENSG[0-9]*");
		Pattern patternEnst = Pattern.compile("ENST[0-9]*");
		Pattern patternRefSeqNm = Pattern.compile("NM_[0-9]*");
		Pattern patternRefSeqNp = Pattern.compile("NP_[0-9]*");

		// Find matching gene IDs for all entities
		int countMatched = 0, countUnMatched = 0;
		for (Entity e : this) {
			String name = e.getName();

			List<String> gids = null;
			String gname = null;

			// Try to match: From most specific to least specific
			if (hasRegexp(patternEnst, name)) {
				// Found ENSEMBL transcript
				gname = findRegexp(patternEnst, name);
				gids = geneIds.getId2tr().getIds(gname);
			} else if (hasRegexp(patternEnsg, name)) {
				// Found ENSEMBL gene ID
				gname = findRegexp(patternEnsg, name);
				gids = new LinkedList<String>();
				gids.add(gname);
			} else if (hasRegexp(patternRefSeqNm, name)) {
				// Found RefeSeq mRNA ID
				gname = findRegexp(patternRefSeqNm, name);
				gids = geneIds.getId2refseqId().getIds(gname);
			} else if (hasRegexp(patternRefSeqNp, name)) {
				// Found RefeSeq protein ID
				gname = findRegexp(patternRefSeqNp, name);
				gids = geneIds.getId2refseqProtId().getIds(gname);
			} else {
				// May be it's gene name followed by some other stuff
				gname = name.split(" ")[0];
				gids = geneIds.getId2geneName().getIds(gname);
			}

			// Found any GeneIDs
			if (gids != null) {
				countMatched++;

				// Add all gene IDs
				for (String gid : gids)
					add(e, gid);

			} else countUnMatched++;

		}

		Timer.showStdErr("Done. Entities matched to geneIDs:" + countMatched + " / " + countUnMatched);
	}

	/**
	 * Load a two-column file into a Hash
	 * @param name
	 * @param fileName
	 * @param map
	 */
	protected void loadMap(String name, String fileName, HashMap<String, String> map) {
		Timer.showStdErr("Loading " + name + " from '" + fileName + "'");

		int i = 1;
		for (String line : Gpr.readFile(fileName).split("\n")) {
			String rec[] = line.split("\t");
			String id = rec[0];
			String componentId = rec[1];
			map.put(id, componentId);
			Gpr.showMark(i++, SHOW_EVERY);
		}

		System.err.println("");
		Timer.showStdErr("Total objects loaded: " + map.size());
	}

	/**
	 * Load pathway events
	 * @param name
	 * @param fileName
	 * @param map
	 */
	protected void loadPathway2HasEvent() {
		String name = "Pathway_2_hasEvent";
		String fileName = dirName + name + ".txt";
		Timer.showStdErr("Loading " + name + " from '" + fileName + "'");

		int i = 1;
		for (String line : Gpr.readFile(fileName).split("\n")) {
			// Parse line
			String rec[] = line.split("\t");
			String id = rec[0];
			String eventId = rec[1];

			if (id.equals("DB_ID")) continue; // Skip title

			// Add event to pathway
			Pathway pathway = (Pathway) getEntity(id);
			Event event = (Event) getEntity(eventId);
			pathway.add(event);

			Gpr.showMark(i++, SHOW_EVERY);
		}

		System.err.println("");
		Timer.showStdErr("Total events assigned: " + (i - 1));
	}

	/**
	 * Load compartment information
	 * @param name
	 * @param fileName
	 * @param map
	 */
	protected void loadPhysicalEntity2Compartment() {
		String name = "PhysicalEntity_2_compartment";
		String fileName = dirName + name + ".txt";
		Timer.showStdErr("Loading " + name + " from '" + fileName + "'");

		int i = 1;
		for (String line : Gpr.readFile(fileName).split("\n")) {
			// Parse line
			String rec[] = line.split("\t");
			String id = rec[0];
			String compartmentId = rec[1];

			if (id.equals("DB_ID")) continue; // Skip title

			// Get entity & compartment
			Entity e = getEntity(id);
			Compartment compartment = (Compartment) getEntity(compartmentId);

			// Assign compartment (if not already assigned)
			if (e.getCompartment() != null) throw new RuntimeException("Compartment already assigned for entity: " + e);
			e.setCompartment(compartment);

			Gpr.showMark(i++, SHOW_EVERY);
		}

		System.err.println("");
		Timer.showStdErr("Total compartments assigned: " + (i - 1));
	}

	/**
	 * Load reaction catalyst
	 * @param name
	 * @param fileName
	 * @param map
	 */
	protected void loadReactionlikeEvent2CatalystActivity() {
		String name = "ReactionlikeEvent_2_catalystActivity";
		String fileName = dirName + name + ".txt";
		Timer.showStdErr("Loading " + name + " from '" + fileName + "'");

		int i = 1;
		for (String line : Gpr.readFile(fileName).split("\n")) {
			// Parse line
			String rec[] = line.split("\t");
			String id = rec[0];
			String catalystId = rec[1];

			if (id.equals("DB_ID")) continue; // Skip title

			// Add event to pathway
			Reaction reaction = (Reaction) entityById.get(id);
			if (reaction == null) continue; // Reaction not found? Skip

			Entity e = getEntity(catalystId);
			reaction.addCatalyst(e);

			Gpr.showMark(i++, SHOW_EVERY);
		}

		System.err.println("");
		Timer.showStdErr("Total outputs assigned: " + (i - 1));
	}

	/**
	 * Load reaction inputs
	 * @param name
	 * @param fileName
	 * @param map
	 */
	protected void loadReactionlikeEvent2Input() {
		String name = "ReactionlikeEvent_2_input";
		String fileName = dirName + name + ".txt";
		Timer.showStdErr("Loading " + name + " from '" + fileName + "'");

		int i = 1;
		for (String line : Gpr.readFile(fileName).split("\n")) {
			// Parse line
			String rec[] = line.split("\t");
			String id = rec[0];
			String inputId = rec[1];

			if (id.equals("DB_ID")) continue; // Skip title

			// Add event to pathway
			Reaction reaction = (Reaction) entityById.get(id);
			if (reaction == null) continue; // Reaction not found? Skip

			Entity e = getEntity(inputId);
			reaction.addInput(e);

			Gpr.showMark(i++, SHOW_EVERY);
		}

		System.err.println("");
		Timer.showStdErr("Total inputs assigned: " + (i - 1));
	}

	/**
	 * Load reaction outputs
	 * @param name
	 * @param fileName
	 * @param map
	 */
	protected void loadReactionlikeEvent2Output() {
		String name = "ReactionlikeEvent_2_output";
		String fileName = dirName + name + ".txt";
		Timer.showStdErr("Loading " + name + " from '" + fileName + "'");

		int i = 1;
		for (String line : Gpr.readFile(fileName).split("\n")) {
			// Parse line
			String rec[] = line.split("\t");
			String id = rec[0];
			String outputId = rec[1];

			if (id.equals("DB_ID")) continue; // Skip title

			// Add event to pathway
			Reaction reaction = (Reaction) entityById.get(id);
			if (reaction == null) continue; // Reaction not found? Skip

			Entity e = getEntity(outputId);
			if (reaction.getId() == 74711) //
				Gpr.debug(reaction.getName() + "\t--->\t" + e.toStringSimple());
			reaction.addOutput(e);

			Gpr.showMark(i++, SHOW_EVERY);
		}

		System.err.println("");
		Timer.showStdErr("Total outputs assigned: " + (i - 1));
	}

	/**
	 * Load reaction regulation
	 * @param name
	 * @param fileName
	 * @param map
	 */
	protected void loadRegulation() {
		String name = "Regulation";
		String fileName = dirName + name + ".txt";
		Timer.showStdErr("Loading " + name + " from '" + fileName + "'");

		int i = 1;
		for (String line : Gpr.readFile(fileName).split("\n")) {
			// Parse line
			String rec[] = line.split("\t");
			String id = rec[0];
			String regulatedEntityId = rec[1];
			String regulatorId = rec[2];

			if (id.equals("DB_ID")) continue; // Skip title

			// Add event to pathway
			//			Gpr.debug("Adding:\tregulatedEntityId: " + regulatedEntityId + "\t" + objectType.get(regulatedEntityId));
			Reaction reaction = (Reaction) entityById.get(regulatedEntityId);
			if (reaction == null) continue; // Reaction not found? Skip

			Entity e = getEntity(regulatorId);
			reaction.addRegulator(e, objectType.get(id));

			Gpr.showMark(i++, SHOW_EVERY);
		}

		System.err.println("");
		Timer.showStdErr("Total regulations assigned: " + (i - 1));
	}

	/**
	 * Reset all nodes in the circuit
	 */
	public void reset() {
		for (Entity e : this)
			e.reset();
	}

	/**
	 * Scale weights
	 */
	void scaleWeights() {
		for (Entity e : this)
			if (e.isReaction()) ((Reaction) e).scaleWeights();
	}

	/**
	 * Set input nodes (fixed outputs from GTEx values)
	 * @param gtex
	 */
	void setInputsFromGtex(GtexExperiment gtexExperiment) {
		Gtex gtex = gtexExperiment.getGtex();

		for (String gid : gtex.getGeneIds()) {
			List<Entity> entities = entitiesByGeneId.get(gid);

			if (entities != null) {
				double value = gtexExperiment.getValue(gid);
				if (!Double.isNaN(value)) {
					for (Entity e : entities)
						e.setFixedOutput(value);
				}
			}
		}
	}

	@Override
	public String toString() {
		CountByType countByType = new CountByType();

		for (Entity e : entityById.values())
			countByType.inc(e.getClass().getSimpleName());

		return countByType.toString();
	}

	/**
	 * Show details
	 * @return
	 */
	public String toStringDetails() {
		StringBuilder sb = new StringBuilder();

		for (Entity e : entityById.values())
			sb.append(e + "\n");

		return sb.toString();
	}

	/**
	 * Run some simulations
	 * @param gtex
	 * @param gtexExperiment
	 */
	public boolean zzz(GtexExperiment gtexExperiment) {
		// Initialize 
		if (monitor == null) monitor = createMonitor(); // Create monitor if needed
		if (monitorTrace == null) monitorTrace = createMonitor("74695"); // Create monitor if needed
		reset(); // Reset previous values
		setInputsFromGtex(gtexExperiment); // Set input nodes (fixed outputs from GTEx values)
		scaleWeights(); // Scale weights

		// Calculate circuit
		calc(gtexExperiment);

		// Add results to monitors
		String experimentLabel = gtexExperiment.getTissueTypeDetail();
		monitor.addResults(experimentLabel);
		monitorTrace.addResults(experimentLabel);

		// Save results
		String file = Gpr.HOME + "/zzz." + monitorTrace.sizeResults() + ".txt";
		Timer.showStdErr("Saving results to '" + file + "'");
		monitor.save(file);

		return true;
	}
}
