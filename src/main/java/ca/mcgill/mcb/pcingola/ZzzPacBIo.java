package ca.mcgill.mcb.pcingola;

import java.util.ArrayList;

import ca.mcgill.mcb.pcingola.fastq.Fastq;
import ca.mcgill.mcb.pcingola.fileIterator.FastqFileIterator;
import ca.mcgill.mcb.pcingola.fileIterator.NeedlemanWunsch;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;

public class ZzzPacBIo {

	public static void main(String[] args) {
		String fastq = Gpr.HOME + "/scyth/pacbio/SCYTH-8Kb_0.38nM_bind-1-2_G02_1.fastq";

		// Read data
		Timer.showStdErr("Reading  file: " + fastq);
		ArrayList<Fastq> seqs = new ArrayList<Fastq>();
		FastqFileIterator fqfi = new FastqFileIterator(fastq);
		for (Fastq fq : fqfi)
			seqs.add(fq);
		Timer.showStdErr("done. Sequences: " + seqs.size());

		// Compare all to all
		int count = 0;
		for (int i = 0; i < seqs.size(); i++) {
			Fastq fq1 = seqs.get(i);
			for (int j = i + 1; j < seqs.size(); j++) {
				Timer.showStdErr("Align :" + i + " , " + j + "\tCount: " + count);
				Fastq fq2 = seqs.get(j);

				NeedlemanWunsch nw = new NeedlemanWunsch(fq1.getSequence(), fq2.getSequence());
				nw.align();
				count++;
			}
		}
	}
}
