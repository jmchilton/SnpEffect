package ca.mcgill.mcb.pcingola.spliceSites;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;

import ca.mcgill.mcb.pcingola.fileIterator.LineFileIterator;
import ca.mcgill.mcb.pcingola.stats.CountByType;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.GprSeq;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * Analyze sequences from splice sites
 * 
 * @author pcingola
 */
public class SpliceSequenceAnalysis {

	public static final double THRESHOLD_ENTROPY = 0.05;
	public static final int THRESHOLD_COUNT = 100;
	public static final double THRESHOLD_P = 0.95;
	public static final int BRANCH_SIZE = 5;
	public static final int ACCEPTOR_SIZE = 5;
	//	public static final int NMER_SIZE = 6;
	//	public static final int NMER_NULL_ITERATIONS = 100 * 1000 * 1000;
	public static final double MIN_OE_BRANCH = 5.0;

	ArrayList<String> donors = new ArrayList<String>();
	ArrayList<String> accs = new ArrayList<String>();
	ArrayList<String> branches = new ArrayList<String>();
	AcgtTree acgtTreeDonors = new AcgtTree();
	AcgtTree acgtTreeAcc = new AcgtTree();
	AcgtTree acgtTreeBranch = new AcgtTree();
	HashMap<String, Integer> donorAcc = new HashMap<String, Integer>();

	double thresholdPDonor;
	double thresholdEntropyDonor;
	double thresholdPAcc;
	double thresholdEntropyAcc;

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

	/**
	 * Find acceptors for this donor
	 * @param donorSeq
	 */
	void acc4donor(String donorSeq) {
		// Create a new tree using all these sequences
		AcgtTree tree = new AcgtTree();
		for (int i = 0; i < donors.size(); i++) {
			String donor = donors.get(i);
			if (donor.startsWith(donorSeq)) {
				String acc = GprSeq.reverse(accs.get(i));
				tree.add(acc);
			}
		}

		// Show them
		for (String accSeq : tree.findNodeNames(thresholdEntropyAcc, thresholdPAcc, THRESHOLD_COUNT)) {
			if (accSeq.length() > 1) {
				accSeq = GprSeq.reverse(accSeq);
				add(donorSeq, accSeq);
			}
		}
	}

	void add(String donor, String acceptor) {
		String key = String.format("%-10s\t%10s", donor, acceptor);
		int count = countDonorAcc(donor, acceptor);
		donorAcc.put(key, count);
	}

	/**
	 * Count how many entries that have both 'donor' and 'acceptor' 
	 * @param donor
	 * @param acceptor
	 * @return
	 */
	int countDonorAcc(String donor, String acceptor) {
		int count = 0;
		for (int i = 0; i < donors.size(); i++) {
			String d = donors.get(i);
			String a = accs.get(i);

			if (d.startsWith(donor) && a.endsWith(acceptor)) count++;
		}
		return count;
	}

	void countMotifMatch(CountByType countMotif, String branchStr) {
		for (String motif : SpliceBranchAnalysis.BRANCH_POINT_SEQUENCES) {
			int idx = branchStr.indexOf(motif);
			if (idx >= 0) countMotif.inc(motif);
		}
	}

	/**
	 * Count branch motifs in entries that have both 'donor' and 'acceptor' 
	 * @param donor
	 * @param acceptor
	 * @return
	 */
	void countMotifMatch(String donor, String acceptor) {
		StringBuilder sb = new StringBuilder();
		StringBuilder fasta = new StringBuilder();

		CountByType countMotif = new CountByType();
		int countBranch = 0, fastaId = 0;
		for (int i = 0; i < donors.size(); i++) {
			String d = donors.get(i);
			String a = accs.get(i);

			if (d.startsWith(donor) && a.endsWith(acceptor)) {
				String branch = branches.get(i);
				countMotifMatch(countMotif, branch);
				countBranch += branch.length();

				fasta.append(">id_" + fastaId + "\n" + branch.subSequence(0, branch.length() - acceptor.length()) + "\n");
				fastaId++;
			}
		}

		boolean show = false;
		for (String m : countMotif.getTypeList()) {
			double expected = countBranch * Math.pow(0.25, m.length());
			double oe = countMotif.get(m) / expected;
			if (oe >= MIN_OE_BRANCH) {
				sb.append(String.format("\t\t\t\t\t\t%-10s\t%d\t%4.2f\n", m, countMotif.get(m), oe));
				show = true;
			}
		}
		if (show) System.out.print("\t\t\t\t\t\tBranch:\n" + sb);

		// Write fasta file 
		String fastaFile = SpliceBranchAnalysis.OUTPUT_DIR + "/" + donor + "-" + acceptor + ".fa";
		Timer.showStdErr("Writing fasta sequences to file: " + fastaFile);
		Gpr.toFile(fastaFile, fasta);
	}

	/**
	 * Find donors for this acceptor
	 * @param accSeq
	 */
	void donor4acc(String accSeq) {
		// Create a new tree using all these sequences
		AcgtTree tree = new AcgtTree();
		for (int i = 0; i < accs.size(); i++) {
			String acc = GprSeq.reverse(accs.get(i));
			if (acc.endsWith(accSeq)) {
				String donor = donors.get(i);
				tree.add(donor);
			}
		}

		// Show them
		for (String donorSeq : tree.findNodeNames(thresholdEntropyDonor, thresholdPDonor, THRESHOLD_COUNT))
			if (donorSeq.length() > 1) add(donorSeq, accSeq);

	}

	/**
	 * Find an probability threshold using THRESHOLD_P quantile
	 * @param tree
	 * @return
	 */
	double findEntropyThreshold(AcgtTree tree) {
		List<Double> values = tree.entropyAll(THRESHOLD_COUNT);
		Collections.sort(values);
		int index = (int) (values.size() * THRESHOLD_ENTROPY);
		return values.get(index);
	}

	/**
	 * Find an probability threshold using THRESHOLD_P quantile
	 * @param tree
	 * @return
	 */
	double findPthreshold(AcgtTree tree) {
		List<Double> values = tree.pAll(THRESHOLD_COUNT);
		Collections.sort(values);
		int index = (int) (values.size() * THRESHOLD_P);
		return values.get(index);
	}

	public void run(String fileName) {
		//---
		// Process file
		//---
		Timer.showStdErr("Reading sequences from " + fileName);
		LineFileIterator lfi = new LineFileIterator(fileName);
		for (String line : lfi) {
			String seqs[] = line.split("\t");

			String donor = seqs[0].toUpperCase();
			String branch = seqs[1].toUpperCase();
			String acc = seqs[2].toUpperCase();

			if ((donor.indexOf('N') < 0) && (branch.indexOf('N') < 0) && (acc.indexOf('N') < 0)) {
				donors.add(donor);
				acgtTreeDonors.add(seqs[0]);

				branches.add(branch);

				accs.add(acc);
				acgtTreeAcc.add(GprSeq.reverse(seqs[2]));
			}
		}
		Timer.showStdErr("Done. Total added: " + donors.size());

		//---
		// Find donor - acceptor pairs
		//---
		Timer.showStdErr("Calculate thresholds");
		thresholdPDonor = findPthreshold(acgtTreeDonors);
		thresholdEntropyDonor = findEntropyThreshold(acgtTreeDonors);
		thresholdPAcc = findPthreshold(acgtTreeAcc);
		thresholdEntropyAcc = findEntropyThreshold(acgtTreeAcc);

		Timer.showStdErr("Donors Thresholds:\tEntropy: " + thresholdEntropyDonor + "\t\tProbability: " + thresholdPDonor);
		for (String seq : acgtTreeDonors.findNodeNames(thresholdEntropyDonor, thresholdPDonor, THRESHOLD_COUNT)) {
			if (seq.length() > 1) acc4donor(seq);
		}

		Timer.showStdErr("Find acceptors");
		Timer.showStdErr("Acceptors Thresholds:\tEntropy: " + thresholdEntropyAcc + "\t\tProbability: " + thresholdPAcc);
		for (String seq : acgtTreeAcc.findNodeNames(thresholdEntropyAcc, thresholdPAcc, THRESHOLD_COUNT)) {
			if (seq.length() > 1) donor4acc(GprSeq.reverse(seq));
		}

		//---
		// Show all donor - acc pairs (sort by number of matches)
		//---
		Timer.showStdErr("Add Donor - Acceptors pairs: ");
		ArrayList<String> keys = new ArrayList<String>();
		keys.addAll(donorAcc.keySet());
		Collections.sort(keys, new Comparator<String>() {

			@Override
			public int compare(String arg0, String arg1) {
				return donorAcc.get(arg1) - donorAcc.get(arg0);
			}
		});

		for (String key : keys) {
			String k[] = key.trim().split("\\s+");
			System.out.println(donorAcc.get(key) + "\t" + key);
			countMotifMatch(k[0], k[1]);
		}

		Timer.showStdErr("Finished!");
	}

}
