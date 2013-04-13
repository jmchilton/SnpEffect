package ca.mcgill.mcb.pcingola;

import ca.mcgill.mcb.pcingola.motif.Jaspar;
import ca.mcgill.mcb.pcingola.motif.Pwm;
import ca.mcgill.mcb.pcingola.util.Gpr;

public class Zzz {

	Jaspar jaspar;
	StringBuilder outSb = new StringBuilder();

	/**
	 * Main
	 * @param args
	 */
	public static void main(String[] args) {
		String dir = Gpr.HOME + "/snpEff/db/jaspar";
		String jasparMatrixFile = dir + "/jaspar.bin";
		String outFile = dir + "/stats.txt";

		Zzz zzz = new Zzz(jasparMatrixFile);
		zzz.run();
		Gpr.toFile(outFile, zzz.outSb.toString());
	}

	public Zzz(String jasparMatrixFile) {
		jaspar = new Jaspar();
		jaspar.load(jasparMatrixFile);
	}

	String changeSequence(char bestSequence[], int pos, char newBase) {
		char newSeq[] = new char[bestSequence.length];

		for (int i = 0; i < bestSequence.length; i++)
			newSeq[i] = bestSequence[i];

		newSeq[pos] = newBase;
		return new String(newSeq);
	}

	void run() {
		// Title
		String out = "pwm\tbestSequence\tbestScore\tpos\tnewBase\tnewSeq\tisConserved\tdiff";
		System.out.println(out);
		outSb.append(out + "\n");

		for (Pwm pwm : jaspar) {
			char bestSeq[] = pwm.bestSequence();
			String bestSequence = new String(bestSeq);
			double bestScore = pwm.score(bestSequence);

			for (int pos = 0; pos < pwm.length(); pos++) {
				for (char newBase : Pwm.BASES) {
					String newSeq = changeSequence(bestSeq, pos, newBase);

					if (!newSeq.equals(bestSequence)) {
						double newScore = pwm.score(newSeq);
						double diff = newScore - bestScore;

						out = pwm.getId() + "\t" + bestSequence + "\t" + bestScore + "\t" + pos + "\t" + newBase + "\t" + newSeq + "\t" + Boolean.toString(pwm.isConserved(pos)).toUpperCase() + "\t" + diff;
						System.out.println(out);
						outSb.append(out + "\n");
					}
				}
			}
		}
	}
}
