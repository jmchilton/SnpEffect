package ca.mcgill.mcb.pcingola.interval;

import ca.mcgill.mcb.pcingola.binseq.DnaSequence;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * Interval for the whole chromosome
 * If a SNP has no 'ChromosomeInterval' => it is outside the chromosome => Invalid
 * 
 * @author pcingola
 *
 */
public class Chromosome extends Marker {

	private static final long serialVersionUID = 1636197649250882952L;

	double chromosomeNum;
	DnaSequence sequence = null;

	/**
	 * Simplify chromosome name
	 * @param chNameOri
	 * @return
	 */
	public static String simpleName(String chrName) {
		return ChromosomeSimpleName.get(chrName);
	}

	public Chromosome(Marker parent, int start, int end, int strand, String id) {
		super(null, start, end, strand, id); // Parent = null to avoid sanity check (it will always fail for chromosomes)
		this.parent = parent;
		type = EffectType.CHROMOSOME.toString();
		setChromosomeName(id);
	}

	/**
	 * Compare only chromosome's name 
	 * @param i2
	 * @return
	 */
	public int compareChromoName(Interval interval) {
		Chromosome i2 = (Chromosome) interval;

		// Compare chromosome
		if ((chromosomeNum == 0) || (i2.chromosomeNum == 0)) return id.compareTo(i2.id); // Use string comparisson

		// Use numeric comparisson
		if (chromosomeNum - i2.chromosomeNum < 0) return -1;
		if (chromosomeNum - i2.chromosomeNum > 0) return 1;
		return 0;
	}

	public DnaSequence getDnaSequence() {
		return sequence;
	}

	public String getSequence() {
		return sequence.toString();
	}

	@Override
	protected boolean isShowWarningIfParentDoesNotInclude() {
		return false;
	}

	/**
	 * Set chromosome name
	 * Note: Removes prefixes (such as 'chr') and parse numeric version.
	 * E.g. 'chr2' becomes '2'. Also numeric '2' is assigned to 'chromosomeNum' to facilitate order by number (so that '2' is ordered before '21')
	 * @param chromo
	 */
	private void setChromosomeName(String chromo) {
		id = simpleName(chromo);
		chromosomeNum = Gpr.parseIntSafe(id); // Try to parse a numeric string
	}

	public void setLength(int len) {
		end = len - 1; // Remember that intervals are zero-based
	}

	/**
	 * Set sequence for this chromosome
	 * @param sequenceStr
	 */
	public void setSequence(String sequenceStr) {
		sequence = new DnaSequence(sequenceStr, true);
		setLength(sequenceStr.length()); // Update chromosome length

	}

}
