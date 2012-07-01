package ca.mcgill.mcb.pcingola.spliceSites;

import java.util.ArrayList;

import ca.mcgill.mcb.pcingola.fileIterator.LineFileIterator;
import ca.mcgill.mcb.pcingola.util.GprSeq;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * Analyze sequences from splice sites
 * 
 * @author pcingola
 */
public class SpliceSequenceAnalysis {

	public static final double THRESHOLD_ENTROPY = 1.5;
	public static final int THRESHOLD_COUNT = 50;
	public static final int BRANCH_SIZE = 5;

	ArrayList<String> donors = new ArrayList<String>();
	ArrayList<String> accs = new ArrayList<String>();
	ArrayList<String> branches = new ArrayList<String>();
	AcgtTree acgtTreeDonors = new AcgtTree("");
	AcgtTree acgtTreeAcc = new AcgtTree("");
	AcgtTree acgtTreeBranch = new AcgtTree("");

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

			String donor = seqs[0].toUpperCase();
			String branch = seqs[1].toUpperCase();
			String acc = seqs[2].toUpperCase();

			if (donor.indexOf('N') < 0) {
				donors.add(donor);
				acgtTreeDonors.add(seqs[0]);
			}

			if (branch.indexOf('N') < 0) {
				branches.add(branch);
				for (int i = 0; i < branch.length() - BRANCH_SIZE; i++) {
					String b = branch.substring(i, i + BRANCH_SIZE + 1);
					acgtTreeBranch.add(b);
				}
			}

			if (acc.indexOf('N') < 0) {
				accs.add(acc);
				acgtTreeAcc.add(GprSeq.reverse(seqs[2]));
			}
		}

		Timer.showStdErr("Done. Total added: " + donors.size());
		System.out.println("\n\nDonors: " + acgtTreeDonors.seqConservation() + "\n" + acgtTreeDonors.toString("", THRESHOLD_ENTROPY, THRESHOLD_COUNT));
		for (String seq : acgtTreeDonors.findNodeNames(THRESHOLD_ENTROPY, THRESHOLD_COUNT))
			System.out.println("\t" + seq);;

		System.out.println("\n\nAcceptors: " + acgtTreeAcc.seqConservation() + "\n" + acgtTreeAcc.toString("", THRESHOLD_ENTROPY, THRESHOLD_COUNT));
		for (String seq : acgtTreeAcc.findNodeNames(THRESHOLD_ENTROPY, THRESHOLD_COUNT))
			System.out.println("\t" + seq);;

		//System.out.println("\n\nBranch:\n" + acgtTreeBranch.toString("", 1.95, THRESHOLD_COUNT));
	}
}
