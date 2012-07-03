package ca.mcgill.mcb.pcingola.spliceSites;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import ca.mcgill.mcb.pcingola.fileIterator.FastaFileIterator;
import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.motif.MotifLogo;
import ca.mcgill.mcb.pcingola.motif.Pwm;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.stats.CountByType;
import ca.mcgill.mcb.pcingola.stats.IntStats;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.GprSeq;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * Analyze sequences from splice sites
 * 
 * @author pcingola
 */
public class SpliceAnalysis {

	/**
	 * A set of PWMs
	 * @author pablocingolani
	 *
	 */
	class PwmSet {
		String name;
		Pwm pwmAcc, pwmDonor;
		CountByType countMotif;
		IntStats lenStats;
		int motifMatchedBases = 0, motifMatchedStr = 0;
		int updates = 0;

		public PwmSet(String name) {
			this.name = name;
			pwmAcc = new Pwm(2 * SIZE_SPLICE + 1);
			pwmDonor = new Pwm(2 * SIZE_SPLICE + 1);
			lenStats = new IntStats();
			countMotif = new CountByType();
		}

		void len(int len) {
			lenStats.sample(len);
		}

		@Override
		public String toString() {
			StringBuilder out = new StringBuilder();

			out.append("<tr>\n");
			out.append("\t<td> <b>" + name + "</b> </td>\n");
			out.append("\t<td> " + updates + "</td>\n");

			// Donor motif
			MotifLogo mlDonor = new MotifLogo(pwmDonor);
			out.append("\t<td>\n");
			out.append(mlDonor.toStringHtml(HTML_WIDTH, HTML_HEIGHT));
			out.append("\t</td>\n");

			// Branch motif match count
			out.append("\t<td>\n\t<table>\n");
			out.append("\t\t<tr> <th> Motif </th> <th> Observed </th> <th> Expected </th> <th> Obs. / Exp. </th> </tr> \n");
			for (String motif : countMotif.getTypeList()) {
				int totalBases = motifMatchedBases - motifMatchedStr * motif.length();
				double expected = totalBases * Math.pow(0.25, motif.length());
				double oe = countMotif.get(motif) / expected; // ratio = Observed / expected

				// Colors
				String bg = "ffffff";
				if (oe > 10) bg = "ff0000";
				else if (oe > 2) bg = "ff8888";
				else if (oe > 1.2) bg = "ffcccc";

				out.append(String.format("\t<tr bgcolor=%s> <td>%s </td> <td align=right> %d </td> <td align=right> %.1f </td> <td align=right> %.2f </td> </tr> \n", bg, motif, countMotif.get(motif), expected, oe));
			}
			out.append("\t</table></center> </td>\n");

			// Acceptor motif
			MotifLogo mlAcc = new MotifLogo(pwmAcc);
			out.append("\t<td>\n");
			out.append(mlAcc.toStringHtml(HTML_WIDTH, HTML_HEIGHT));
			out.append("\t</td>\n");

			out.append("</tr>\n");

			return out.toString();
		}

		public void update(String accStr, String donorStr) {
			updates++;
			if (accStr != null) pwmAcc.updateCounts(accStr);
			if (donorStr != null) pwmDonor.updateCounts(donorStr);
		}
	}

	public static final String OUTPUT_DIR = Gpr.HOME + "/snpEff/splice";
	public static int SIZE_SPLICE = 10;
	public static int SIZE_CONSENSUS_DONOR = 2;
	public static int SIZE_CONSENSUS_ACCEPTOR = 2;
	public static int SIZE_BRANCH = 60;
	public static final double THRESHOLD_ENTROPY = 0.05;
	public static final int THRESHOLD_COUNT = 100;
	public static final double THRESHOLD_P = 0.95;
	public static final double THRESHOLD_BRANCH_U12_SCORE = 0.99;
	public static final int BRANCH_SIZE = 5;
	public static final int ACCEPTOR_SIZE = 5;
	public static final double MIN_OE_BRANCH = 5.0;
	public static int HTML_WIDTH = 20;
	public static int HTML_HEIGHT = 100;
	public static double MIN_UPDATES_PERC = 0.0000; // Don't show if there are less than this number on the whole genome

	String genomeVer;
	String genomeFasta;
	StringBuilder out = new StringBuilder();

	Config config;

	ArrayList<String> donorsList = new ArrayList<String>();
	ArrayList<String> acceptorsList = new ArrayList<String>();
	ArrayList<String> branchesList = new ArrayList<String>();

	ArrayList<String> daPairDonor = new ArrayList<String>();
	ArrayList<String> daPairAcc = new ArrayList<String>();

	AcgtTree acgtTreeDonors = new AcgtTree();
	AcgtTree acgtTreeAcc = new AcgtTree();

	Pwm pwmU12;
	HashMap<String, Integer> donorAcc = new HashMap<String, Integer>();
	HashMap<String, PwmSet> pwms = new HashMap<String, PwmSet>();

	double thresholdPDonor;
	double thresholdEntropyDonor;
	double thresholdPAcc;
	double thresholdEntropyAcc;
	double thresholdU12Score;

	int countIntrons = 0;

	public static void main(String[] args) {
		// Command line arguments
		if (args.length != 1) {
			System.err.println("Usage: SpliceBranchAnalysis2 genomeName");
			System.exit(1);
		}

		// Find dono-acceptor combinations
		String genomeName = args[0];
		SpliceAnalysis splBr = new SpliceAnalysis(genomeName);
		splBr.run();
	}

	public SpliceAnalysis(String genomeVer) {
		Config config = new Config(genomeVer);
		this.genomeVer = genomeVer;

		// Find genome reference
		genomeFasta = config.getFileNameGenomeFasta();
		if (genomeFasta == null) throw new RuntimeException("Cannot find reference genome: " + config.getFileListGenomeFasta());
	}

	/**
	 * Find acceptors for this donor
	 * @param donorSeq
	 */
	void acc4donor(String donorSeq) {
		// Create a new tree using all these sequences
		AcgtTree tree = new AcgtTree();
		for (int i = 0; i < donorsList.size(); i++) {
			String donor = donorsList.get(i);
			if (donor.startsWith(donorSeq)) {
				String acc = GprSeq.reverse(acceptorsList.get(i));
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
	 * Find the best score for PWM matrix in U12 branch points
	 * @param seq
	 * @return
	 */
	double bestU12Score(String seq) {
		int max = seq.length() - pwmU12.length();
		double best = 0;
		for (int i = 0; i < max; i++) {
			String sub = seq.substring(i, i + pwmU12.length());
			double score = pwmU12.score(sub);
			best = Math.max(best, score);
		}

		return best;
	}

	/**
	 * Calculate threshold of U12 PWM scores  
	 */
	void branchU12Threshold() {
		ArrayList<Double> scores = new ArrayList<Double>();

		for (String branch : branchesList) {
			double bestScore = bestU12Score(branch);
			scores.add(bestScore);
		}

		// Get quantile
		Collections.sort(scores);
		int index = (int) (THRESHOLD_BRANCH_U12_SCORE * scores.size());
		thresholdU12Score = scores.get(index);
	}

	/**
	 * Count how many entries that have both 'donor' and 'acceptor' 
	 * @param donor
	 * @param acceptor
	 * @return
	 */
	int countDonorAcc(String donor, String acceptor) {
		int count = 0;
		for (int i = 0; i < donorsList.size(); i++) {
			String d = donorsList.get(i);
			String a = acceptorsList.get(i);

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
	 * Find donors for this acceptor
	 * @param accSeq
	 */
	void donor4acc(String accSeq) {
		// Create a new tree using all these sequences
		AcgtTree tree = new AcgtTree();
		for (int i = 0; i < acceptorsList.size(); i++) {
			String acc = GprSeq.reverse(acceptorsList.get(i));
			if (acc.endsWith(accSeq)) {
				String donor = donorsList.get(i);
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

	PwmSet getPwmSet(String key) {
		PwmSet ps = pwms.get(key);
		if (ps == null) {
			Gpr.debug("Creating PWMs for key: '" + key + "'");
			ps = new PwmSet(key);
			pwms.put(key, ps);
		}
		return ps;
	}

	/**
	 * Lad data from files
	 */
	void load() {
		Timer.showStdErr("Loading U12 PWM: " + genomeVer);
		pwmU12 = new Pwm(OUTPUT_DIR + "/u12_branch.txt");

		Timer.showStdErr("Loading: " + genomeVer);
		config = new Config(genomeVer);
		config.loadSnpEffectPredictor();
		Timer.showStdErr("done");
	}

	/**
	 * Show and append an output line
	 * @param line
	 */
	void out(Object o) {
		String s = o.toString();
		out.append(s + "\n");
		//System.out.println(s);
	}

	public void run() {
		// Load data
		load();

		// Find splice sequences
		spliceSequences();

		//		// Find donor acceptor pairs
		//		spliceDonoAcceptorPairs();
		//		acgtTreeDonors = acgtTreeAcc = null; // Free some unused objects

		branchU12Threshold();

		//		// Splice site PWM analysis
		//		spliceAnalysis();
		//
		//		//---
		//		// Save output
		//		//---
		//		String outputFile = OUTPUT_DIR + "/" + this.getClass().getSimpleName() + "_" + genomeVer + ".html";
		//		Timer.showStdErr("Saving output to: " + outputFile);
		//		Gpr.toFile(outputFile, out);
	}

	void spliceAnalysis() {
		out("<pre>\n");

		// Iterate over all chromosomes
		FastaFileIterator ffi = new FastaFileIterator(genomeFasta);
		for (String chrSeq : ffi) {
			String chrName = Chromosome.simpleName(ffi.getName());
			Timer.showStdErr("Reading chromosome: " + chrName);
			spliceAnalysis(chrName, chrSeq);
		}
		out("</pre>\n");

		// Show PwmSets
		int threshold = (int) (MIN_UPDATES_PERC * countIntrons);
		Gpr.debug("Exons: " + countIntrons + "\tThreshold: " + threshold);
		ArrayList<String> names = new ArrayList<String>();
		names.addAll(pwms.keySet());
		Collections.sort(names);
		out("<table border=1>\n");
		out("<tr> <th> Donor type </th>  <th> Count </th>  <th> Donor Motif </th> <th> Branch matches </th> <th> Acceptor Motif </th> <th> Branch (best energy) Motif </th> <th> Intron length </th> </tr>\n");
		for (String key : names)
			if (pwms.get(key).updates > threshold) out(pwms.get(key));
		out("</table>\n");
	}

	/**
	 * Run analysis for one chromosome
	 * @param chrName
	 * @param chrSeq
	 */
	void spliceAnalysis(String chrName, String chrSeq) {
		int countEx = 0;
		HashSet<String> done = new HashSet<String>();

		for (Gene gene : config.getGenome().getGenes()) {
			if (gene.getChromosomeName().equals(chrName)) { // Same chromosome
				for (Transcript tr : gene) {
					int prev = -1;
					for (Exon ex : tr.sorted()) {
						countEx++;

						if (prev >= 0) {
							if (prev > ex.getStart()) System.err.println("WARNINIG: Exon check failed. Skipping: " + ex);
							else {
								int start = prev;
								int end = ex.getStart();
								String key = chrName + ":" + start + "-" + end;

								if (!done.contains(key)) updatePwm(tr, chrSeq, start, end);

								done.add(key);
							}
						}

						prev = ex.getEnd();
					}
				}
			}
		}

		Timer.showStdErr("Chromo: " + chrName + "\tGenes: " + config.getGenome().getGenes().size() + "\tExons: " + countEx);
	}

	void spliceDonoAcceptorPairs() {
		//---
		// Create trees
		//---
		Timer.showStdErr("Create quaternary trees");
		for (String donor : donorsList)
			if (donor.indexOf('N') < 0) acgtTreeDonors.add(donor);

		for (String acc : acceptorsList)
			if (acc.indexOf('N') < 0) acgtTreeAcc.add(GprSeq.reverse(acc));

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
			String da[] = key.trim().split("\\s+");
			daPairDonor.add(da[0]);
			daPairAcc.add(da[1]);

			System.out.println(donorAcc.get(key) + "\t" + key);
		}

		Timer.showStdErr("Finished!");
	}

	/**
	 * Find splice sequences for this genome
	 */
	void spliceSequences() {
		// Iterate over all chromosomes
		FastaFileIterator ffi = new FastaFileIterator(genomeFasta);
		for (String chrSeq : ffi) {
			String chrName = Chromosome.simpleName(ffi.getName());
			Timer.showStdErr("Reading chromosome: " + chrName);
			spliceSequences(chrName, chrSeq);
		}
	}

	/**
	 * Find splice sequences for this cromosome
	 * @param chrName
	 * @param chrSeq
	 */
	void spliceSequences(String chrName, String chrSeq) {
		int countEx = 0;
		HashSet<String> done = new HashSet<String>();

		for (Gene gene : config.getGenome().getGenes()) {
			if (gene.getChromosomeName().equals(chrName)) { // Same chromosome
				for (Transcript tr : gene) {
					int prev = -1;
					for (Exon ex : tr.sorted()) {
						countEx++;

						if (prev >= 0) {
							if (prev > ex.getStart()) System.err.println("WARNINIG: Exon check failed. Skipping: " + ex);
							else {
								int start = prev;
								int end = ex.getStart();
								String key = chrName + ":" + start + "-" + end;

								// Already added? (do not add twice)
								if (!done.contains(key)) spliceSequences(tr, chrSeq, start, end);

								done.add(key);
							}
						}

						prev = ex.getEnd();
					}
				}
			}
		}

		Timer.showStdErr("Chromo: " + chrName + "\tGenes: " + config.getGenome().getGenes().size() + "\tExons: " + countEx);
	}

	/**
	 * Find splice sequences for this intron
	 */
	void spliceSequences(Transcript tr, String chrSeq, int intronStart, int intronEnd) {
		countIntrons++;

		//---
		// Get donor and acceptor coordinates and strings
		//---
		int splDonorStart = intronStart - SIZE_SPLICE;
		int splDonorEnd = intronStart + SIZE_SPLICE;
		int splAccStart = intronEnd - SIZE_SPLICE;
		int splAccEnd = intronEnd + SIZE_SPLICE;
		int splBranchStart = intronEnd - SIZE_BRANCH + 1;
		int splBranchEnd = intronEnd;

		// Get strings
		if (tr.isStrandMinus()) {
			splDonorStart = intronEnd - SIZE_SPLICE;
			splDonorEnd = intronEnd + SIZE_SPLICE;
			splAccStart = intronStart - SIZE_SPLICE;
			splAccEnd = intronStart + SIZE_SPLICE;
			splBranchStart = intronStart + 1;
			splBranchEnd = intronStart + SIZE_BRANCH;
		}

		String donorStr = chrSeq.substring(splDonorStart, splDonorEnd + 1).toUpperCase();
		String accStr = chrSeq.substring(splAccStart, splAccEnd + 1).toUpperCase();
		String branchStr = chrSeq.substring(splBranchStart, splBranchEnd).toUpperCase();

		if (tr.isStrandMinus()) {
			donorStr = GprSeq.reverseWc(donorStr);
			accStr = GprSeq.reverseWc(accStr);
			branchStr = GprSeq.reverseWc(branchStr);
		}

		String intronSeqDonor = donorStr.substring(SIZE_SPLICE + 1);
		String intronSeqAcc = accStr.substring(0, SIZE_SPLICE);

		// Add to arrays
		donorsList.add(intronSeqDonor);
		acceptorsList.add(intronSeqAcc);
		branchesList.add(branchStr);
	}

	/**
	 * Update PWM
	 * @param tr
	 * @param intronStart
	 * @param intronEnd
	 */
	void updatePwm(Transcript tr, String chrSeq, int intronStart, int intronEnd) {
		countIntrons++;

		int len = intronEnd - intronStart;

		//---
		// Get donor and acceptor coordinates and strings
		//---
		int splDonorStart = intronStart - SIZE_SPLICE;
		int splDonorEnd = intronStart + SIZE_SPLICE;
		int splAccStart = intronEnd - SIZE_SPLICE;
		int splAccEnd = intronEnd + SIZE_SPLICE;
		int splBranchStart = intronEnd - SIZE_BRANCH + 1;
		int splBranchEnd = intronEnd;

		// Get strings
		if (tr.isStrandMinus()) {
			splDonorStart = intronEnd - SIZE_SPLICE;
			splDonorEnd = intronEnd + SIZE_SPLICE;
			splAccStart = intronStart - SIZE_SPLICE;
			splAccEnd = intronStart + SIZE_SPLICE;
			splBranchStart = intronStart + 1;
			splBranchEnd = intronStart + SIZE_BRANCH;
		}

		String donorStr = chrSeq.substring(splDonorStart, splDonorEnd + 1).toUpperCase();
		String accStr = chrSeq.substring(splAccStart, splAccEnd + 1).toUpperCase();
		String branchStr = chrSeq.substring(splBranchStart, splBranchEnd).toUpperCase();

		if (tr.isStrandMinus()) {
			donorStr = GprSeq.reverseWc(donorStr);
			accStr = GprSeq.reverseWc(accStr);
			branchStr = GprSeq.reverseWc(branchStr);
		}

		String intronSeqDonor = donorStr.substring(SIZE_SPLICE + 1);
		String intronSeqAcc = accStr.substring(0, SIZE_SPLICE);

		//---
		// Group by donor type
		//---

		// Donor consensus ('GT' or 'GC'?)
		String donorConsensus = donorStr.substring(SIZE_SPLICE + 1, SIZE_SPLICE + 1 + SIZE_CONSENSUS_DONOR);
		if (donorConsensus.indexOf('N') >= 0) return; // Ignore if there is an 'N'

		// Use long consensus? U12
		String accConsensus = accStr.substring(SIZE_SPLICE - SIZE_CONSENSUS_ACCEPTOR, SIZE_SPLICE);
		if (donorConsensus.indexOf('N') >= 0) return; // Ignore if there is an 'N'
		String consensus = donorConsensus + "_" + accConsensus;

		int maxLenDa = 0;
		for (int i = 0; i < daPairDonor.size(); i++) {
			String don = daPairDonor.get(i);
			String ac = daPairAcc.get(i);
			if (intronSeqDonor.startsWith(don) && intronSeqAcc.endsWith(ac)) {
				int lenda = don.length() + ac.length();
				if (lenda > maxLenDa) {
					maxLenDa = lenda;
					donorConsensus = don;
					accConsensus = ac;
				}
			}
		}

		//---
		// Update PWM
		//---
		PwmSet pwmSet = getPwmSet(consensus);
		pwmSet.update(accStr, donorStr);
		pwmSet.len(len);

		// Update total counts
		pwmSet = getPwmSet(" ALL");
		pwmSet.update(accStr, donorStr);
		pwmSet.len(len);
	}

}
