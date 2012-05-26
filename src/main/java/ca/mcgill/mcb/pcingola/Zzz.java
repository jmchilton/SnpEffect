package ca.mcgill.mcb.pcingola;

import ca.mcgill.mcb.pcingola.fileIterator.FastaFileIterator;
import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.motif.MotifLogo;
import ca.mcgill.mcb.pcingola.motif.Pwm;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.stats.CountByType;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.GprSeq;
import ca.mcgill.mcb.pcingola.util.Timer;

public class Zzz {

	public static int SIZE_SPLICE = 2;
	public static int SIZE_BRANCH = 60;
	public static int POLYPYRIMIDINE_TRACT_SIZE = 3;
	public static boolean debug = false;

	public static int HTML_WIDTH = 20;
	public static int HTML_HEIGHT = 100;

	public static String[] MOTIFS = { "TCCTTAAC" // U12 consensus, from "Evolitionary Fates and Origins of U12 
			, "TACTAAC" // From Yeast
			, "CTAAC", "CTAAT", "CTCAC", "CTGAC", "CTGAT" // CTRAY Motif
	};

	int motifMatchedBases = 0, motifMatchedStr = 0;
	String genomeVer;
	String genomeFasta;
	Config config;
	Pwm pwmAcc, pwmDonor, pwmBranch;
	StringBuilder out = new StringBuilder();
	CountByType countMotif = new CountByType();

	public static void main(String[] args) {
		Zzz zzz = new Zzz("testHg3763Chr20");
		//		Zzz zzz = new Zzz("hg19");
		zzz.run();
	}

	public Zzz(String genomeVer) {
		this.genomeVer = genomeVer;
		genomeFasta = Gpr.HOME + "/snpEff/data/genomes/" + genomeVer + ".fa.gz";
		pwmAcc = new Pwm(2 * SIZE_SPLICE + 1);
		pwmDonor = new Pwm(2 * SIZE_SPLICE + 1);
		pwmBranch = new Pwm(2 * POLYPYRIMIDINE_TRACT_SIZE + 1);

		load();
	}

	/**
	 * Calculate and show score
	 * @param branchStr
	 */
	String countMotifMatch(String branchStr) {
		for (String motif : MOTIFS) {
			int idx = branchStr.indexOf(motif);
			if (idx >= 0) countMotif.inc(motif);
		}

		motifMatchedBases += branchStr.length();
		motifMatchedStr++;

		return ""; // String.format("%.1f\t%s\n", maxScore, sb.toString());
	}

	void load() {
		Timer.showStdErr("Loading");
		config = new Config(genomeVer);
		config.loadSnpEffectPredictor();
		Timer.showStdErr("done");
	}

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

		// Motif counts
		out.append("Count Motifs:<br>\n<pre>\n");
		out.append(countMotif + "\n");
		out.append("Expected:\n");
		for (String motif : MOTIFS) {
			double expected = (motifMatchedBases - (motif.length() * motifMatchedStr)) * Math.pow(0.25, motif.length());
			out.append(motif + "\t" + expected + "\n");
		}
		out.append("</pre>\n");

		// Show PWMs
		MotifLogo mlDonor = new MotifLogo(pwmDonor);
		MotifLogo mlAcc = new MotifLogo(pwmAcc);
		MotifLogo mlBranch = new MotifLogo(pwmBranch);

		out.append("Donor:<p>\n");
		out.append(mlDonor.toStringHtml(HTML_WIDTH, HTML_HEIGHT));
		out.append("Acceptor:<p>\n");
		out.append(mlAcc.toStringHtml(HTML_WIDTH, HTML_HEIGHT));
		out.append("Branch:<p>\n");
		out.append(mlBranch.toStringHtml(HTML_WIDTH, HTML_HEIGHT));

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
		// Get donor and acceptor coordinates
		int splDonorStart = intronStart - SIZE_SPLICE;
		int splDonorEnd = intronStart + SIZE_SPLICE;
		int splAccStart = intronEnd - SIZE_SPLICE;
		int splAccEnd = intronEnd + SIZE_SPLICE;
		int splBranchStart = intronEnd - SIZE_BRANCH + 1;
		int splBranchEnd = intronEnd - 2;

		if (tr.isStrandMinus()) {
			splDonorStart = intronEnd - SIZE_SPLICE;
			splDonorEnd = intronEnd + SIZE_SPLICE;
			splAccStart = intronStart - SIZE_SPLICE;
			splAccEnd = intronStart + SIZE_SPLICE;
			splBranchStart = intronStart + 2;
			splBranchEnd = intronStart + SIZE_BRANCH - 1;
		}

		String donorStr = chrSeq.substring(splDonorStart, splDonorEnd + 1);
		String accStr = chrSeq.substring(splAccStart, splAccEnd + 1);
		String branchStr = chrSeq.substring(splBranchStart, splBranchEnd + 1);

		if (tr.isStrandMinus()) {
			donorStr = GprSeq.reverseWc(donorStr);
			accStr = GprSeq.reverseWc(accStr);
			branchStr = GprSeq.reverseWc(branchStr);
		}

		if (debug) {
			String line = (tr.isStrandPlus() ? "+" : "-") + "\tIntron: [" + intronStart + ", " + intronEnd + "]\tlen: " + (intronEnd - intronStart) + "\tDonor: " + donorStr + "\tAcc:" + accStr + "\tBranch: " + branchStr;
			System.out.println(line);
			out.append(line + "\n");
		}

		String scoreStr = countMotifMatch(branchStr);
		out.append(scoreStr);

		// Update PWM
		pwmAcc.updateCounts(accStr);
		pwmDonor.updateCounts(donorStr);
		//pwmBranch.updateCounts(branchStr);
	}
}
