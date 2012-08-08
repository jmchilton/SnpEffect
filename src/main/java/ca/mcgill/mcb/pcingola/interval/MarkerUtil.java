package ca.mcgill.mcb.pcingola.interval;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import ca.mcgill.mcb.pcingola.fileIterator.LineFileIterator;

/**
 * Generic utility methods for Markers
 * 
 * @author pcingola
 */
public class MarkerUtil {

	/**
	 * Read intervals from a file using a simplt TXT format
	 * Format: 
	 * 		chr \t start \t end \t id
	 * 
	 * Note: Zero-based positions
	 * 
	 * @param fileName : Path to file
	 * @param genome : Genome to use. Can be null (a new one will be created)
	 * @param positionBase : Position offset. Use '1' for one-based coordinates and '0' for zero-based coordinates.
	 */
	public static Markers readTxt(String fileName, Genome genome, int positionBase) {
		if (genome == null) genome = new Genome();
		Markers markers = new Markers();

		// Parse lines
		LineFileIterator lfi = new LineFileIterator(fileName);
		int lineNum = 1;
		for (String line : lfi) {
			Marker interval = new Marker();
			interval.readTxt(line, lineNum, genome, positionBase);
			markers.add(interval);
			lineNum++;
		}
		return markers;
	}

	/**
	 * Redundant markers in a list: Find intervals that are totally included in other intervals in the list
	 * @param markersOri
	 * @return A map  markerIncluded -> markerLarge, where  markerIncluded in completely included in markerLarge
	 */
	public static Map<Marker, Marker> redundant(Collection<? extends Marker> markersOri) {
		Map<Marker, Marker> redundant = new HashMap<Marker, Marker>();

		// Find which markers are redundant?
		ArrayList<Marker> markers = new ArrayList<Marker>();
		markers.addAll(markersOri);
		int size = markers.size();

		// Iterate on all markers
		for (int i = 0; i < size; i++) {
			Marker mi = markers.get(i);

			// Is marker 'mi' included in any other marker?
			Marker markerLarge = null;
			for (int j = 0; (j < size) && (markerLarge == null); j++) {
				Marker mj = markers.get(j);
				if ((i != j) && (mj.includes(mi))) { // Not the same interval and it is fully included? 
					if (mi.includes(mj) && (i > j)) {
						// If they are included both ways, it means that they are exactly the same.
						// We have to avoid deleting both of them twice, so we arbitrarely don't add them if (i > j) 
					} else markerLarge = mj;
				}
			}

			// Add to redundant marker
			if (markerLarge != null) redundant.put(mi, markerLarge);
		}

		return redundant;
	}

}
