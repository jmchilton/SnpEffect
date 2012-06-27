package ca.mcgill.mcb.pcingola;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import ca.mcgill.mcb.pcingola.fileIterator.FastaFileIterator;
import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Markers;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.interval.Utr3prime;
import ca.mcgill.mcb.pcingola.interval.Utr5prime;
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
		Pwm pwmAcc, pwmDonor, pwmBestEnergy;
		CountByType countMotif;
		IntStats lenStats;
		int motifMatchedBases = 0, motifMatchedStr = 0;
		int updates = 0;
		IntStats dgStatsAll, dgStatsBest;

		public PwmSet(String name) {
			this.name = name;
			pwmAcc = new Pwm(2 * SIZE_SPLICE + 1);
			pwmDonor = new Pwm(2 * SIZE_SPLICE + 1);
			pwmBestEnergy = new Pwm(SIZE_SPLICE + 1);
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
			int bestDg = 0, bestIdx = 0;
			double bestPpt = 0;
			String bestDgSeq = "";
			int max = branchStr.length() - intron.length();
			for (int i = 0; i < max; i++) {
				String branchSub = branchStr.substring(i, i + intron.length());
				int dg = rnaFreeEnergy.energyInt(intron, branchSub);
				dgStatsAll.sample(dg); // Stats for all Delta_G

				// Lowest energy?
				if (dg <= bestDg) {
					double ppt = polyPyrimidineScore(branchStr, i + intron.length());

					// Better PPT score?
					if (ppt > bestPpt) {
						bestDg = dg;
						bestDgSeq = branchSub;
						bestIdx = branchStr.length() - i;
						bestPpt = ppt;
					}
				}
			}

			// Filter what we show
			if ((bestDg < energyThreshold) && (bestPpt > 0.7)) {
				Gpr.debug("Energy :  " + bestDg + "\tIntronSeq: " + intron + "\tBestSeq: " + bestDgSeq + "\tPPT score: " + bestPpt + "\tIdx: " + bestIdx);
				pwmBestEnergy.updateCounts(bestDgSeq);
			}
			dgStatsBest.sample(bestDg); // Stats for best Delta_G

			return bestDg;
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

		/**
		 * Calculate Number of Polypyrimidines in a region
		 * @param bases
		 * @param start
		 * @param end
		 * @return
		 */
		double polyPyrimidineScore(String branchStr, int start) {
			int score = 0, count = 0;

			if (start < 0) start = 0;
			int end = start + POLYPYRIMIDINE_TRACT_SIZE;

			char bases[] = branchStr.toCharArray();
			for (int i = start; (i <= end) && (i < bases.length); i++) {
				switch (bases[i]) {
				case 'C':
				case 'T':
					score++;
					break;
				}
				count++;
			}

			return score / ((double) count);
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

			// Branch (best energy) motif
			MotifLogo mlBbe = new MotifLogo(pwmBestEnergy);
			out.append("\t<td>\n");
			out.append(mlBbe.toStringHtml(HTML_WIDTH, HTML_HEIGHT));
			out.append("\t</td>\n");

			out.append("\t<td nowrap>\n");
			out.append(dgStatsAll.toString());
			out.append("<br><img src=\"" + dgStatsAll.toStringPlot("Free energy all", "Free Energy", true) + "\">");
			out.append("\t</td>\n");

			out.append("\t<td nowrap>\n");
			out.append(dgStatsBest.toString());
			out.append("<br><img src=\"" + dgStatsBest.toStringPlot("Free energy best", "Free Energy", true) + "\">");
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

	public static final String DELTA_G_FILE = Gpr.HOME + "/workspace/SnpEff/data/stack.deltaG";

	public static boolean debug = false;
	public static boolean test = false;
	public static int SIZE_SPLICE = 9;
	public static int SIZE_BRANCH = 60;
	public static int POLYPYRIMIDINE_TRACT_SIZE = 15;
	public static int BRANCH_ENERGY_LEN = 8;
	public static int MIN_PWM_LEN = 5;
	public static int MAX_PWM_LEN = 10;
	public static double MIN_UPDATES_PERC = 0.0005; // Don't show if there are less than this number on the whole genome
	public static int HTML_WIDTH = 20;
	public static int HTML_HEIGHT = 100;
	public static int ENERGY_DISTRIBUTION_ITERATIONS = 10 * 1000 * 1000;
	public static final double ENERGY_THRESHOLD_QUANTILE = 0.05;

	public static String[] MOTIFS = { "TCCTTAAC", "TCCTTGAC", "TCCTTAAT", "TCCTTGAT" // U12 consensus, from "Evolutionary Fates and Origins of U12 
			, "TACTAAC" // From Yeast
			, "CTAACT", "CTAATT", "CTCACT", "CTGACT", "CTGATT" // 6-letter motifs (CTRAY Motif or CTRAYT?
			, "CTGAC", "CTAAT", "CTGAT", "CTAAC", "CTCAC" // 5-letter motifs
			, "TTCAC" //
	};

	public static String[] DONOR_LONG_KEYS = { "ATATCCT", "GTATCCT" //
			, "GTAAGT", "GTGAGT", "GTAAAA", "GTGAGA" // dbStep
	};

	int countIntrons = 0;
	int energyThreshold = 0;
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
		// Command line arguments
		if (args.length != 1) {
			System.err.println("Usage: SpliceBranchAnalysis genomeName");
			System.exit(1);
		}

		// Run
		String genomeName = args[0];
		SpliceBranchAnalysis splBr = new SpliceBranchAnalysis(genomeName);
		splBr.run();
	}

	public SpliceBranchAnalysis(String genomeVer) {
		Config config = new Config(genomeVer);
		this.genomeVer = genomeVer;

		// Find genome reference
		genomeFasta = config.getFileNameGenomeFasta();
		if (genomeFasta == null) throw new RuntimeException("Cannot find reference genome: " + config.getFileListGenomeFasta());

		pwms = new HashMap<String, PwmSet>();
		load();
	}

	/**
	 * Calculate empirical values for energy and select a threshold
	 * @param quantile
	 * @return
	 */
	int energyThreshold(double quantile) {
		Timer.showStdErr("Energy  statistics");
		IntStats distr = rnaFreeEnergy.empiricalDistribution(BRANCH_ENERGY_LEN, ENERGY_DISTRIBUTION_ITERATIONS);
		int th = (int) distr.getQuantile(quantile);
		Timer.showStdErr("Done. Quantile: " + quantile + "\tThreshold = " + th);
		return th;
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

		rnaFreeEnergy = new RnaStackedFreeEnergy(DELTA_G_FILE);

		Timer.showStdErr("done");
	}

	/**
	 * Show and append an output line
	 * @param line
	 */
	void out(Object o) {
		String s = o.toString();
		out.append(s + "\n");
		System.out.println(s);
	}

	void outStats(IntStats stats, String name) {
		out("<p><center><b>" + name + " sizes</b></center><p>");
		out("<img src=\"" + stats.toStringPlot(name + " sizes", "Size", true) + "\"><br>\n");
		out("<pre>\n");
		out(stats.toString());
		out(stats.toStringHisto());
		out("</pre>\n");
	}

	/**
	 * Run algorithms
	 */
	void run() {
		runGeneLengthAnalysis();
		runSpliceAnalysis();

		//---
		// Save output
		//---
		String outputFile = this.getClass().getSimpleName() + "_" + genomeVer + ".html";
		Timer.showStdErr("Saving output to: " + outputFile);
		Gpr.toFile(outputFile, out);
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
					int prev = -1;
					for (Exon ex : tr.sorted()) {
						countEx++;

						if (prev >= 0) {
							if (prev > ex.getStart()) System.err.println("WARNINIG: Exon check failed. Skipping: " + ex);
							else updatePwm(tr, chrSeq, prev, ex.getStart());
						}

						prev = ex.getEnd();
					}
				}
			}
		}

		Timer.showStdErr("Chromo: " + chrName + "\tGenes: " + config.getGenome().getGenes().size() + "\tExons: " + countEx);
	}

	/**
	 * Perform basic statistics on gene, exon, intron sizes (etc.)
	 */
	void runGeneLengthAnalysis() {
		IntStats sgene = new IntStats();
		IntStats str = new IntStats();
		IntStats sexons = new IntStats();
		IntStats su3 = new IntStats();
		IntStats su5 = new IntStats();
		IntStats sintron = new IntStats();

		for (Gene g : config.getGenome().getGenes()) {
			sgene.sample(g.size());

			for (Transcript tr : g) {
				str.sample(tr.size());

				Markers mtr = new Markers();
				mtr.add(tr);
				Markers mex = new Markers();
				for (Exon ex : tr) {
					sexons.sample(ex.size());
					mex.add(ex);
				}

				for (Utr3prime u3 : tr.get3primeUtrs())
					su3.sample(u3.size());

				for (Utr5prime u5 : tr.get5primeUtrs())
					su5.sample(u5.size());

				// Get introns by substracting exons fomr transcript
				Markers introns = mtr.minus(mex);
				for (Marker m : introns)
					sintron.sample(m.size());
			}
		}

		outStats(sgene, "Gene");
		outStats(str, "Transcript");
		outStats(sexons, "Exon");
		outStats(sintron, "Intron");
		outStats(su3, "UTR 3 prime");
		outStats(su5, "UTR 5 prime");
	}

	void runSpliceAnalysis() {
		// Sample energy distribution
		energyThreshold = energyThreshold(ENERGY_THRESHOLD_QUANTILE / SIZE_BRANCH);

		out("<pre>\n");

		// Iterate over all chromosomes
		FastaFileIterator ffi = new FastaFileIterator(genomeFasta);
		for (String chrSeq : ffi) {
			String chrName = Chromosome.simpleName(ffi.getName());
			Timer.showStdErr("Reading chromosome: " + chrName);
			run(chrName, chrSeq);
		}
		out("</pre>\n");

		// Show PwmSets
		ArrayList<String> names = new ArrayList<String>();
		names.addAll(pwms.keySet());
		Collections.sort(names);
		out("<table border=1>\n");
		out("<tr> <th> Donor type </th>  <th> Count </th>  <th> Donor Motif </th> <th> Branch matches </th> <th> Acceptor Motif </th> <th> Intron length </th> </tr>\n");
		int threshold = (int) (MIN_UPDATES_PERC * countIntrons);
		Gpr.debug("Exons: " + countIntrons + "\tThreshold: " + threshold);
		for (String key : names)
			if (pwms.get(key).updates > threshold) out(pwms.get(key));
		out("</table>\n");
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
			out(line);
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
