package ca.mcgill.mcb.pcingola;

import java.util.HashMap;

import ca.mcgill.mcb.pcingola.fastq.Fastq;
import ca.mcgill.mcb.pcingola.fileIterator.FastqFileIterator;
import ca.mcgill.mcb.pcingola.fileIterator.LineFileIterator;
import ca.mcgill.mcb.pcingola.fileIterator.NeedlemanWunsch;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;

public class ZzzPacBIo {

	String fastq;
	String putativeOverlaps;
	HashMap<String, String> seqByName;

	public static void main(String[] args) {
		String fastq = Gpr.HOME + "/scyth/pacbio/SCYTH-8Kb_0.38nM_bind-1-2_G02_1.fastq";
		String putativeOverlaps = Gpr.HOME + "/scyth/pacbio/SCYTH-8Kb_0.38nM_bind-1-2_G02_1.subreads.100.sort.filter.overlaps.txt";

		ZzzPacBIo zzz = new ZzzPacBIo(fastq, putativeOverlaps);
		zzz.readSeqs();
		zzz.readOverlaps();
	}

	public ZzzPacBIo(String fastq, String putativeOverlaps) {
		this.fastq = fastq;
		this.putativeOverlaps = putativeOverlaps;
		seqByName = new HashMap<String, String>();
	}

	/**
	 * Overlap two sequences
	 * @param fq1
	 * @param fq2
	 */
	void overlap(String readName1, String readName2) {
		Timer.showStdErr("Align :" + readName1 + " , " + readName2);
		String seq1 = seqByName.get(readName1);
		String seq2 = seqByName.get(readName2);
		NeedlemanWunsch nw = new NeedlemanWunsch(seq1, seq2);
		nw.align();
	}

	/**
	 * Read overlap file
	 */
	void readOverlaps() {
		Timer.showStdErr("Reading overlaps file: " + putativeOverlaps);
		LineFileIterator lfi = new LineFileIterator(putativeOverlaps);
		for (String line : lfi) {
			String fields[] = line.split("\t");
			String readName1 = fields[0];

			System.out.println(line);
			System.out.println(readName1);
			for (int i = 1; i < fields.length; i++) {
				String sf[] = fields[i].split(";");
				String readName2 = sf[0]; // Read name

				sf = sf[1].split("/");
				int positiveMaps = Gpr.parseIntSafe(sf[0]); // Positive strand 
				int negativeMaps = sf.length > 1 ? Gpr.parseIntSafe(sf[1]) : 0; // Negative strand

				System.out.println("\t" + readName2 + "\t" + positiveMaps + "\t" + negativeMaps);
				if (positiveMaps > 0) {
					overlap(readName1, readName2);
				}
			}
		}
	}

	/**
	 * Read all sequences
	 */
	void readSeqs() {
		Timer.showStdErr("Reading seuqnece (FASTQ) file: " + fastq);

		FastqFileIterator fqfi = new FastqFileIterator(fastq);
		for (Fastq fq : fqfi) {
			String name = fq.getDescription().substring(1);
			// System.out.println("NAME:" + name);
			seqByName.put(name, fq.getSequence());
		}

		Timer.showStdErr("done. Sequences: " + seqByName.size());
	}
}
