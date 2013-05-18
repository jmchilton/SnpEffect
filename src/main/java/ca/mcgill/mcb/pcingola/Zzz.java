package ca.mcgill.mcb.pcingola;

import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.SnpEffectPredictor;
import ca.mcgill.mcb.pcingola.util.Timer;

public class Zzz {

	public static void main(String[] args) {
		Timer.showStdErr("Loading");

		Config config = new Config("testHg3770Chr22");
		//		Config config = new Config("GRCh37.70");
		SnpEffectPredictor sep = config.loadSnpEffectPredictor();

		Timer.showStdErr("Checking");

		System.out.println(sep.getGenome());

		//		int errorProteinLength = 0;
		//		int errorProteinStopCodons = 0;
		//		int countTranscriptsProteinCoding = 0;
		//
		//		for (Gene g : sep.getGenome().getGenes()) {
		//			for (Transcript tr : g) {
		//
		//				if (tr.isProteinCoding()) countTranscriptsProteinCoding++;
		//
		//				if (tr.isErrorProteinLength()) {
		//					errorProteinLength++; // Protein length error
		//					System.out.println(tr.getId() + "\t" + tr.isProteinCoding() + "\tError: Length\t" + tr.cds().length() + "\t" + tr.protein());
		//				}
		//
		//				if (tr.isErrorStopCodonsInCds()) {
		//					errorProteinStopCodons++; // Protein length error
		//					errorProteinLength++; // Protein length error
		//					System.out.println(tr.getId() + "\t" + tr.isProteinCoding() + "\tError: STOP_IN_CDS\n\t" + tr.cds().length() + "\t" + tr.protein());
		//					System.out.println(tr);
		//				}
		//
		//				if (tr.isErrorStopCodon()) {
		//					System.out.println(tr.getId() + "\t" + tr.isProteinCoding() + "\tError: STOP_CODON\t" + tr.cds().length() + "\t" + tr.protein());
		//				}
		//
		//				if (tr.isErrorStartCodon()) {
		//					System.out.println(tr.getId() + "\t" + tr.isProteinCoding() + "\tError: START_CODON\t" + tr.cds().length() + "\t" + tr.protein());
		//				}
		//
		//			}
		//		}
		//
		//		System.out.println("# Protein coding transcripts : " + countTranscriptsProteinCoding);
		//		System.out.println(String.format("#              Length errors : %06d ( %.2f%% ) ", errorProteinLength, (100.0 * errorProteinLength / countTranscriptsProteinCoding)));
		//		System.out.println(String.format("#          STOP codon errors : %06d ( %.2f%% )", errorProteinStopCodons, (100.0 * errorProteinStopCodons / countTranscriptsProteinCoding)));

		Timer.showStdErr("Done");
	}
}
