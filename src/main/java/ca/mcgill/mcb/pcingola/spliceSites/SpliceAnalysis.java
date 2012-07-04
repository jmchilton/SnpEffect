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
import ca.mcgill.mcb.pcingola.snpEffect.commandLine.SnpEff;
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
public class SpliceAnalysis extends SnpEff {

	/**
	 * A set of PWMs
	 * @author pablocingolani
	 *
	 */
	class PwmSet implements Comparable<PwmSet> {
		String name;
		Pwm pwmAcc, pwmDonor;
		CountByType countMotif;
		IntStats lenStats;
		int motifMatchedBases = 0, motifMatchedStr = 0;
		int updates = 0;
		int countU12 = 0;

		public PwmSet(String name) {
			this.name = name;
			pwmAcc = new Pwm(2 * SIZE_SPLICE + 1);
			pwmDonor = new Pwm(2 * SIZE_SPLICE + 1);
			lenStats = new IntStats();
			countMotif = new CountByType();
		}

		@Override
		public int compareTo(PwmSet ps) {
			int diff = ps.updates - updates;
			if (diff != 0) return diff;
			return name.compareTo(ps.name);
		}

		void incU12() {
			countU12++;
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

			// U12 count
			double expected = updates * (1.0 - THRESHOLD_BRANCH_U12_SCORE);
			double oe = countU12 / expected; // ratio = Observed / expected

			// U12 Colors
			String bg = "ffffff";
			if (oe > 10) bg = "ff0000";
			else if (oe > 2) bg = "ff8888";
			else if (oe > 1.2) bg = "ffcccc";
			out.append(String.format("\t<td bgcolor=%s> <center> %d (%1.2f)" + " </center> </td>\n", bg, countU12, oe));

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

	public static int SIZE_SPLICE = 8;
	public static int SIZE_CONSENSUS_DONOR = 2;
	public static int SIZE_CONSENSUS_ACCEPTOR = 2;
	public static int SIZE_BRANCH = 60;
	public static final double THRESHOLD_ENTROPY = 0.05;
	public static final int THRESHOLD_COUNT = 100;
	public static final double THRESHOLD_P = 0.95;
	public static final double THRESHOLD_BRANCH_U12_SCORE = 0.95;
	public static final int BRANCH_SIZE = 5;
	public static final int ACCEPTOR_SIZE = 5;
	public static final double MIN_OE_BRANCH = 5.0;
	public static int HTML_WIDTH = 20;
	public static int HTML_HEIGHT = 100;
	public static double MIN_UPDATES_PERC = 0.0001; // Don't show if there are less than this number on the whole genome

	String outputDir = ".";
	String genomeVer;
	String genomeFasta;
	StringBuilder out = new StringBuilder();

	Config config;

	ArrayList<String> donorsList = new ArrayList<String>();
	ArrayList<String> acceptorsList = new ArrayList<String>();
	ArrayList<String> branchesList = new ArrayList<String>();

	ArrayList<String> donorAccPairDonor = new ArrayList<String>();
	ArrayList<String> donorAccPairAcc = new ArrayList<String>();

	AcgtTree acgtTreeDonors = new AcgtTree();
	AcgtTree acgtTreeAcc = new AcgtTree();

	Pwm pwmU12;
	HashMap<String, Integer> donorAcc = new HashMap<String, Integer>();
	HashMap<String, PwmSet> pwmSetsByName = new HashMap<String, PwmSet>();

	double thresholdPDonor;
	double thresholdEntropyDonor;
	double thresholdPAcc;
	double thresholdEntropyAcc;
	double thresholdU12Score;

	int countIntrons = 0;

	public SpliceAnalysis() {
		super();
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
				if (acc.indexOf('N') < 0) tree.add(acc);
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
			if (sub.indexOf('N') < 0) {
				double score = pwmU12.score(sub);
				best = Math.max(best, score);
			}
		}

		return best;
	}

	/**
	 * Calculate threshold of U12 PWM scores  
	 */
	void branchU12Threshold() {
		Timer.showStdErr("Finding U12 PWM score distribution and threshold.");
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

	/**
	 * Create one fasta file for each donor-acceptor pair
	 */
	void createSpliceFasta() {
		Timer.showStdErr("Creating FASTA files for each dono-acceptor pair.");

		for (int i = 0; i < donorAccPairDonor.size(); i++) {
			String donor = donorAccPairDonor.get(i);
			String acceptor = donorAccPairAcc.get(i);
			createSpliceFasta(donor, acceptor);
		}
	}

	/**
	 * Count branch motifs in entries that have both 'donor' and 'acceptor' 
	 * @param donor
	 * @param acceptor
	 * @return
	 */
	void createSpliceFasta(String donor, String acceptor) {
		StringBuilder fasta = new StringBuilder();

		int fastaId = 0;
		for (int i = 0; i < donorsList.size(); i++) {
			String d = donorsList.get(i);
			String a = acceptorsList.get(i);

			if (d.startsWith(donor) && a.endsWith(acceptor)) {
				String branch = branchesList.get(i);
				fasta.append(">id_" + fastaId + "\n" + branch.subSequence(0, branch.length() - acceptor.length()) + "\n");
				fastaId++;
			}
		}

		// Write fasta file 
		String fastaFile = outputDir + "/" + genomeVer + "." + donor + "-" + acceptor + ".fa";
		Timer.showStdErr("\tWriting fasta sequences to file: " + fastaFile);
		Gpr.toFile(fastaFile, fasta);
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
				if (donor.indexOf('N') < 0) tree.add(donor);
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
		PwmSet ps = pwmSetsByName.get(key);
		if (ps == null) {
			ps = new PwmSet(key);
			pwmSetsByName.put(key, ps);
		}
		return ps;
	}

	/**
	 * Initialize
	 */
	void init() {
		Timer.showStdErr("Initializing");
		config = new Config(genomeVer);
		genomeFasta = config.getFileNameGenomeFasta();
		if (genomeFasta == null) throw new RuntimeException("Cannot find reference genome: " + config.getFileListGenomeFasta());

		outputDir = config.getDirData() + "/spliceSites";

		// Load data
		load();
	}

	/**
	 * Lad data from files
	 */
	void load() {
		String u12file = config.getDirData() + "/spliceSites/u12_branch.pwm";
		Timer.showStdErr("Loading U12 PWM form file '" + u12file + "'");
		pwmU12 = new Pwm(u12file);

		Timer.showStdErr("Loading: " + genomeVer);
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

	@Override
	public void parseArgs(String[] args) {
		if (args.length != 1) usage("Missing genome");
		genomeVer = args[0];
	}

	@Override
	public boolean run() {
		init();

		// Find splice sequences
		spliceSequences();

		// Find donor acceptor pairs
		spliceDonoAcceptorPairs();
		acgtTreeDonors = acgtTreeAcc = null; // Free some unused objects

		branchU12Threshold();
		createSpliceFasta();

		// Splice site PWM analysis
		spliceAnalysis();

		//---
		// Save output
		//---
		String outputFile = outputDir + "/" + this.getClass().getSimpleName() + "_" + genomeVer + ".html";
		Timer.showStdErr("Saving output to: " + outputFile);
		Gpr.toFile(outputFile, out);

		Timer.showStdErr("Finished!");
		return true;
	}

	void spliceAnalysis() {
		Timer.showStdErr("Splice analysis (PWM). Reading fasta file: " + genomeFasta);

		out("<pre>\n");

		// Iterate over all chromosomes
		FastaFileIterator ffi = new FastaFileIterator(genomeFasta);
		for (String chrSeq : ffi) {
			String chrName = Chromosome.simpleName(ffi.getName());
			spliceAnalysis(chrName, chrSeq);
		}
		out("</pre>\n");

		// Show PwmSets
		int threshold = (int) (MIN_UPDATES_PERC * countIntrons);
		Timer.showStdErr("Filter out low count splice sites. Exons: " + countIntrons + "\tThreshold: " + threshold);
		ArrayList<PwmSet> pwmsets = new ArrayList<PwmSet>();
		pwmsets.addAll(pwmSetsByName.values());
		Collections.sort(pwmsets);
		out("<table border=1>\n");
		out("<tr> <th> Donor type </th>  <th> Count </th>  <th> Donor Motif </th> <th> U12 matches (Observed / Expected) </th> <th> Acceptor Motif </th> </tr>\n");
		for (PwmSet pwmset : pwmsets)
			if (pwmset.updates > threshold) out(pwmset);
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

		Timer.showStdErr("\tChromosome: " + chrName + "\tGenes: " + config.getGenome().getGenes().size() + "\tExons: " + countEx);
	}

	void spliceDonoAcceptorPairs() {
		//---
		// Create trees
		//---
		Timer.showStdErr("Finding donor-acceptor pairs: Creating quaternary trees");
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

		Timer.showStdErr("Donors Thresholds:\n\t\t\tEntropy: " + thresholdEntropyDonor + "\n\t\t\tProbability: " + thresholdPDonor);
		for (String seq : acgtTreeDonors.findNodeNames(thresholdEntropyDonor, thresholdPDonor, THRESHOLD_COUNT)) {
			if (seq.length() > 1) acc4donor(seq);
		}

		Timer.showStdErr("Find acceptors");
		Timer.showStdErr("Acceptors Thresholds:\n\t\t\tEntropy: " + thresholdEntropyAcc + "\n\t\t\tProbability: " + thresholdPAcc);
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
			donorAccPairDonor.add(da[0]);
			donorAccPairAcc.add(da[1]);

			Timer.showStdErr("\t\t" + donorAcc.get(key) + "\t" + key);
		}
	}

	/**
	 * Find splice sequences for this genome
	 */
	void spliceSequences() {
		Timer.showStdErr("Finding splice sequences. Reading fasta file: " + genomeFasta);

		// Iterate over all chromosomes
		FastaFileIterator ffi = new FastaFileIterator(genomeFasta);
		for (String chrSeq : ffi) {
			String chrName = Chromosome.simpleName(ffi.getName());
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

		Timer.showStdErr("\tChromosome: " + chrName + "\tGenes: " + config.getGenome().getGenes().size() + "\tExons: " + countEx);
	}

	/**
	 * Find splice sequences for this intron
	 */
	void spliceSequences(Transcript tr, String chrSeq, int intronStart, int intronEnd) {
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

		int maxLenDa = 0;
		for (int i = 0; i < donorAccPairDonor.size(); i++) {
			String don = donorAccPairDonor.get(i);
			String ac = donorAccPairAcc.get(i);
			if (intronSeqDonor.startsWith(don) && intronSeqAcc.endsWith(ac)) {
				int lenda = don.length() + ac.length();
				if (lenda > maxLenDa) {
					maxLenDa = lenda;
					donorConsensus = don;
					accConsensus = ac;
				}
			}
		}
		String consensus = donorConsensus + "_" + accConsensus;

		//---
		// Branch U12 score
		//---
		double bu12score = bestU12Score(branchStr);

		//---
		// Update PWM
		//---
		PwmSet pwmSet = getPwmSet(consensus);
		pwmSet.update(accStr, donorStr);
		pwmSet.len(len);
		if (bu12score >= thresholdU12Score) pwmSet.incU12();

		// Update total counts
		pwmSet = getPwmSet(" ALL");
		pwmSet.update(accStr, donorStr);
		pwmSet.len(len);

	}

	@Override
	public void usage(String message) {
		if (message != null) System.err.println("Error: " + message + "\n");
		System.err.println("Usage: snpEff genome_version");
		System.exit(-1);
	}

}
