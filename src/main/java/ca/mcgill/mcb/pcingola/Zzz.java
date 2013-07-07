package ca.mcgill.mcb.pcingola;

import java.util.HashSet;
import java.util.List;

import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.interval.Utr;
import ca.mcgill.mcb.pcingola.interval.Utr5prime;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.SnpEffectPredictor;
import ca.mcgill.mcb.pcingola.stats.CountByType;
import ca.mcgill.mcb.pcingola.util.Timer;

public class Zzz {

	public static final int READ_LENGTH = 50;
	public static final int NUM_READS_MARKER = 5;

	Config config;
	SnpEffectPredictor sep;

	public static void main(String[] args) {
		String genome = "testHg3766Chr1";
		genome = "GRCh37.71";
		Zzz zzz = new Zzz(genome);
		zzz.run();
	}

	public Zzz(String genome) {
		config = new Config(genome);
	}

	public void run() {
		Timer.showStdErr("Loading");
		sep = config.loadSnpEffectPredictor();
		Timer.showStdErr("Building");
		sep.buildForest();
		Timer.showStdErr("Done");

		CountByType countByKozak = new CountByType();

		HashSet<String> done = new HashSet<String>();
		for (Gene gene : sep.getGenome().getGenes()) {

			for (Transcript tr : gene) {
				if (!tr.isProteinCoding()) continue;
				if (tr.hasErrorOrWarning()) continue;

				// Get UTR 'key' and check we haven't used this position before
				List<Utr5prime> utrs5 = tr.get5primeUtrs();
				int pos = tr.isStrandPlus() ? 0 : Integer.MAX_VALUE;
				for (Utr utr : utrs5) {
					if (tr.isStrandPlus()) pos = Math.max(pos, utr.getEnd());
					else pos = Math.min(pos, utr.getStart());
				}
				String key = tr.getChromosomeName() + ":" + pos;
				if (done.contains(key)) continue;
				done.add(key);

				// Get UTR sequence
				if (utrs5.size() <= 0) continue;
				Utr5prime utr5 = utrs5.get(0);
				String utr5Str = utr5.getSequence().toLowerCase();
				String cds = tr.cds();

				if ((utr5Str.length() >= 10) && (cds.length() >= 5)) {
					String kozak = utr5Str.substring(utr5Str.length() - 6) + cds.substring(0, 4).toUpperCase();
					System.out.println(kozak + "\t\t\t" + gene.getGeneName() + "\t" + tr.getId() + "\t" + key);
					countByKozak.inc(kozak);
				}
			}
		}

		System.out.println("\t\t\tTop: \n" + countByKozak.toStringTop(100));

	}
}
