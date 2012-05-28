package ca.mcgill.mcb.pcingola;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

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

public class SpliceBranchAnalysis {

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
		IntStats dgStatsAll, dgStatsBest;

		public PwmSet(String name) {
			this.name = name;
			pwmAcc = new Pwm(2 * SIZE_SPLICE + 1);
			pwmDonor = new Pwm(2 * SIZE_SPLICE + 1);
			lenStats = new IntStats();
			countMotif = new CountByType();

			dgStatsAll = new IntStats();
			dgStatsBest = new IntStats();

		}

		/**
		 * Calculate the best free energy potition
		 * 
		 * @param intron
		 * @param branchStr
		 * @return
		 */
		public double bestFreeEergy(String intron, String branchStr) {
			if ((intron.indexOf('N') >= 0) || (branchStr.indexOf('N') >= 0)) return 0;

			// RWC
			intron = GprSeq.reverseWc(intron).toUpperCase();

			// Energy on each subsequence of the branch site
			double bestDg = 0;
			String bestBranch = "";
			int max = branchStr.length() - intron.length();
			for (int i = 0; i < max; i++) {
				String branchSub = branchStr.substring(i, i + intron.length());
				double dg = rnaFreeEnergy.energy(intron, branchSub);
				dgStatsAll.sample((int) (dg * 10)); // Stats for all Delta_G

				// Lowest energy?
				if (dg < bestDg) {
					bestDg = dg;
					bestBranch = branchSub;
				}
			}

			dgStatsBest.sample((int) (bestDg * 10)); // Stats for best Delta_G
			return 0.0;
		}

		/**
		 * Calculate and show score
		 * @param branchStr
		 */
		void countMotifMatch(String branchStr) {
			for (String motif : MOTIFS) {
				int idx = branchStr.indexOf(motif);
				if (idx >= 0) countMotif.inc(motif);
			}

			motifMatchedBases += branchStr.length();
			motifMatchedStr++;
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

			out.append("\t<td>\n");
			out.append(lenStats.toString());
			out.append("<img src=\"" + dgStatsAll.toStringPlot("Free energy all", "Free Energy", true) + "\">");
			out.append("\t</td>\n");

			out.append("\t<td>\n");
			out.append(lenStats.toString());
			out.append("<img src=\"" + dgStatsBest.toStringPlot("Free energy best", "Free Energy", true) + "\">");
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

	public static boolean debug = false;
	public static boolean test = false;
	public static int SIZE_SPLICE = 9;
	public static int SIZE_BRANCH = 60;
	public static int POLYPYRIMIDINE_TRACT_SIZE = 3;
	public static int BRANCH_ENERGY_LEN = 8;
	public static int MIN_PWM_LEN = 5;
	public static int MAX_PWM_LEN = 10;
	public static double MIN_UPDATES_PERC = 0.0005; // Don't show if there are less than this number on the whole genome
	public static int HTML_WIDTH = 20;
	public static int HTML_HEIGHT = 100;

	public static String[] MOTIFS = { "TCCTTAAC", "TCCTTGAC", "TCCTTAAT", "TCCTTGAT" // U12 consensus, from "Evolutionary Fates and Origins of U12 
			, "TACTAAC" // From Yeast
			, "CTAACT", "CTAATT", "CTCACT", "CTGACT", "CTGATT" // CTRAY Motif or CTRAYT?
	};

	public static String[] DONOR_LONG_KEYS = { "ATATCCT", "GTATCCT" };

	int countIntrons = 0;
	String genomeVer;
	String genomeFasta;
	StringBuilder out = new StringBuilder();
	Config config;
	HashMap<String, PwmSet> pwms;
	RnaStackedFreeEnergy rnaFreeEnergy;

	/**
	 * Main
	 * @param args
	 */
	public static void main(String[] args) {

		SpliceBranchAnalysis zzz = new SpliceBranchAnalysis(test ? "testHg3763Chr20" // 
				: "hg19");
		//		: "GRCh37.66");

		zzz.run();
	}

	public SpliceBranchAnalysis(String genomeVer) {
		this.genomeVer = genomeVer;
		genomeFasta = Gpr.HOME + "/snpEff/data/genomes/" + genomeVer + ".fa.gz";
		pwms = new HashMap<String, PwmSet>();
		load();
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
		Timer.showStdErr("Loading: " + genomeVer);
		config = new Config(genomeVer);
		config.loadSnpEffectPredictor();

		rnaFreeEnergy = new RnaStackedFreeEnergy(Gpr.HOME + "/workspace/SnpEff/data/stack.deltaG");

		Timer.showStdErr("done");
	}

	/**
	 * Run algorithms
	 */
	void run() {
		out.append("<pre>\n");

		// Iterate over all chromosomes
		FastaFileIterator ffi = new FastaFileIterator(genomeFasta);
		for (String chrSeq : ffi) {
			String chrName = Chromosome.simpleName(ffi.getName());
			Timer.showStdErr("Reading chromosome: " + chrName);
			run(chrName, chrSeq);
		}
		out.append("</pre>\n");

		// Show PwmSets
		ArrayList<String> names = new ArrayList<String>();
		names.addAll(pwms.keySet());
		Collections.sort(names);
		out.append("<table border=1>\n");
		//out.append("<tr> <th> Donor type </th>  <th> Count </th>  <th> Donor Motif </th> <th> Branch matches </th> <th> Branch K-mers </th>  <th> Acceptor Motif </th> </tr>\n");
		out.append("<tr> <th> Donor type </th>  <th> Count </th>  <th> Donor Motif </th> <th> Branch matches </th> <th> Acceptor Motif </th> <th> Intron length </th> </tr>\n");
		int threshold = (int) (MIN_UPDATES_PERC * countIntrons);
		Gpr.debug("Exons: " + countIntrons + "\tThreshold: " + threshold);
		for (String key : names)
			if (pwms.get(key).updates > threshold) out.append(pwms.get(key));
		out.append("</table>\n");

		System.out.println(out);
		Gpr.toFile(Gpr.HOME + "/zzz.html", out);
	}

	/**
	 * Run analysis for one chromosome
	 * @param chrName
	 * @param chrSeq
	 */
	void run(String chrName, String chrSeq) {
		int countEx = 0;

		for (Gene gene : config.getGenome().getGenes()) {
			if (gene.getChromosomeName().equals(chrName)) { // Same chromosome
				for (Transcript tr : gene) {
					if (tr.isProteinCoding()) { // Is protein coding?
						int prev = -1;
						for (Exon ex : tr.sorted()) {
							countEx++;

							if (prev >= 0) {
								if (prev > ex.getStart()) throw new RuntimeException("WTF!");
								updatePwm(tr, chrSeq, prev, ex.getStart());
							}

							prev = ex.getEnd();
						}
					}
				}
			}
		}

		Timer.showStdErr("Chromo: " + chrName + "\tGenes: " + config.getGenome().getGenes().size() + "\tExons: " + countEx);
	}

	/**
	 * Calculate Number of Polypyrimidines in a region
	 * @param bases
	 * @param start
	 * @param end
	 * @return
	 */
	double score(char bases[], int start, int end) {
		int score = 0, count = 0;
		if (start < 0) start = 0;
		for (int i = start; (i <= end) && (i < bases.length); i++) {
			switch (bases[i]) {
			case 'A':
			case 'G':
				score--;
				break;

			case 'C':
			case 'T':
				score++;
				break;

			default:
				throw new RuntimeException("Unknown base '" + bases[i] + "'");
			}
			count++;
		}

		if (count < POLYPYRIMIDINE_TRACT_SIZE) return 0;
		return score / ((double) count);
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
			splBranchStart = intronStart;
			splBranchEnd = intronStart + SIZE_BRANCH - 1;
		}

		String donorStr = chrSeq.substring(splDonorStart, splDonorEnd + 1).toUpperCase();
		String accStr = chrSeq.substring(splAccStart, splAccEnd + 1).toUpperCase();
		String branchStr = chrSeq.substring(splBranchStart, splBranchEnd + 1).toUpperCase();

		if (tr.isStrandMinus()) {
			donorStr = GprSeq.reverseWc(donorStr);
			accStr = GprSeq.reverseWc(accStr);
			branchStr = GprSeq.reverseWc(branchStr);
		}

		//---
		// Group by donor type
		//---

		// Donor consensus ('GT' or 'GC'?)
		String donorConsensus = donorStr.substring(SIZE_SPLICE + 1, SIZE_SPLICE + 3);
		if (donorConsensus.indexOf('N') >= 0) return; // Ignore if there is an 'N'
		// Use long consensus? U12
		String donorConsensusLong = donorStr.substring(SIZE_SPLICE + 1);
		for (String dlk : DONOR_LONG_KEYS)
			if (donorConsensusLong.startsWith(dlk)) donorConsensus = dlk;

		String intronSeq = donorStr.substring(SIZE_SPLICE + 1);

		if (debug) {
			String line = (tr.isStrandPlus() ? "+" : "-") + "\tIntron: [" + intronStart + ", " + intronEnd + "]\tlen: " + (intronEnd - intronStart) + "\tDonor: " + donorStr + "\tAcc:" + accStr + "\tBranch: " + branchStr;
			System.out.println(line);
			out.append(line + "\n");
		}

		//---
		// Update PWM
		//---
		PwmSet pwmSet = getPwmSet(donorConsensus);
		pwmSet.update(accStr, donorStr);
		pwmSet.countMotifMatch(branchStr);
		pwmSet.bestFreeEergy(intronSeq, branchStr);
		pwmSet.len(len);

		// Update total counts
		pwmSet = getPwmSet(" ALL");
		pwmSet.update(accStr, donorStr);
		pwmSet.countMotifMatch(branchStr);
		pwmSet.bestFreeEergy(intronSeq, branchStr);
		pwmSet.len(len);
	}
}
