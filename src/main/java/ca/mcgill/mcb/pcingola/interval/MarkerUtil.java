package ca.mcgill.mcb.pcingola.interval;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;

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
	 * @return A list of markers that are totally included in other markers in the list
	 */
	public static Collection<Marker> redundant(Collection<? extends Marker> markersOri) {
		HashSet<Marker> redundant = new HashSet<Marker>();

		// Find which markers are redundant?
		ArrayList<Marker> markers = new ArrayList<Marker>();
		markers.addAll(markersOri);
		int size = markers.size();

		// Iterate on all markers
		for (int i = 0; i < size; i++) {
			Marker mi = markers.get(i);

			// Is marker 'mi' included in any other marker?
			boolean included = false;
			for (int j = 0; (j < size) && (!included); j++) {
				Marker mj = markers.get(j);
				if ((i != j) && (!redundant.contains(mj))) {
					included = mj.includes(mi);
				}
			}

			if (included) redundant.add(mi); // Add to redundant list
		}

		return redundant;
	}

}
