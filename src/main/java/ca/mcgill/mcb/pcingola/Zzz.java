package ca.mcgill.mcb.pcingola;

import java.util.List;

import ca.mcgill.mcb.pcingola.fileIterator.SeqChangeBedFileIterator;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * Simple test program
 * @author pcingola
 */
public class Zzz {

	public static void main(String[] args) {
		String bedFile = "/Users/pablocingolani/askat/bed/GRCh37.69.genes.bed.gz";

		Timer.showStdErr("Loading: " + bedFile);
		SeqChangeBedFileIterator bed = new SeqChangeBedFileIterator(bedFile);
		List<SeqChange> intervals = bed.load();

		Timer.showStdErr("Done: " + intervals.size());
		for (Marker m : intervals) {
			String mid = m.getId().replaceAll("[^a-zA-Z0-9\\-\\.]+", "_");
			Gpr.debug(m + "\t" + m.getId() + "\t" + mid);
		}
	}
}
