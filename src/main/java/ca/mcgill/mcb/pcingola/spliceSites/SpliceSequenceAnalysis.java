package ca.mcgill.mcb.pcingola.spliceSites;

import java.util.ArrayList;

import ca.mcgill.mcb.pcingola.fileIterator.LineFileIterator;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * Analyze sequences from splice sites
 * 
 * @author pcingola
 */
public class SpliceSequenceAnalysis {

	ArrayList<String> donors = new ArrayList<String>();
	ArrayList<String> accs = new ArrayList<String>();
	ArrayList<String> branches = new ArrayList<String>();
	AcgtTree acgtTreeDonors = new AcgtTree("");

	public static void main(String[] args) {
		// Command line arguments
		if (args.length != 1) {
			System.err.println("Usage: SpliceBranchAnalysis2 spliceFile");
			System.exit(1);
		}

		// Run
		String fileName = args[0];
		SpliceSequenceAnalysis splBr = new SpliceSequenceAnalysis();
		splBr.run(fileName);
	}

	public void run(String fileName) {
		Timer.showStdErr("Reading sequences from " + fileName);
		LineFileIterator lfi = new LineFileIterator(fileName);
		for (String line : lfi) {
			String seqs[] = line.split("\t");

			donors.add(seqs[0]);
			branches.add(seqs[1]);
			accs.add(seqs[2]);

			acgtTreeDonors.add(seqs[0]);
		}

		Timer.showStdErr("Done. Total added: " + donors.size());
		System.out.println(acgtTreeDonors.toString("", 120));
	}
}
