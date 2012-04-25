package ca.mcgill.mcb.pcingola.interval;

import java.util.List;

import ca.mcgill.mcb.pcingola.codons.CodonTable;
import ca.mcgill.mcb.pcingola.codons.CodonTables;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * An interval intended as a mark
 * 
 * @author pcingola
 */
public class Marker extends Interval {

	private static final long serialVersionUID = 7878886900660027549L;
	protected String type = "";

	public Marker(Marker parent, int start, int end, int strand, String id) {
		super(parent, start, end, strand, id);

		// Adjust parent if child is not included?
		if ((parent != null) && !parent.includes(this)) {
			String err = "";
			if (isShowWarningIfParentDoesNotInclude()) err = "WARNING: " + this.getClass().getSimpleName() + " is not included in parent " + parent.getClass().getSimpleName() + ". " //
					+ "\t" + this.getClass().getSimpleName() + " '" + this.getId() + "'  [ " + this.getStart() + " , " + this.getEnd() + " ]" //
					+ "\t" + parent.getClass().getSimpleName() + " '" + parent.getId() + "' [ " + parent.getStart() + " , " + parent.getEnd() + " ]";

			// Adjust parent?
			if (isAdjustIfParentDoesNotInclude(parent)) {
				parent.adjust(this);
				if (isShowWarningIfParentDoesNotInclude()) err += "\t=> Adjusting " + parent.getClass().getSimpleName() + " '" + parent.getId() + "' to [ " + parent.getStart() + " , " + parent.getEnd() + " ]";
			}

			// Show an error?
			if (isShowWarningIfParentDoesNotInclude()) System.err.println(err);
		}
	}

	/**
	 * Adjust [start,end] to include child
	 * @param child
	 */
	protected void adjust(Marker child) {
		start = Math.min(start, child.getStart());
		end = Math.max(end, child.getEnd());
	}

	/**
	 * Get a suitable codon table
	 * @return
	 */
	public CodonTable codonTable() {
		return CodonTables.getInstance().getTable(getGenome(), getChromosomeName());
	}

	/**
	 * Compare by start and end
	 */
	@Override
	public int compareTo(Interval i2) {
		// Compare chromosome names
		Marker m2 = (Marker) i2;

		Chromosome chr1 = getChromosome();
		if (chr1 != null) {
			Chromosome chr2 = m2.getChromosome();
			if (chr2 != null) {
				int compChromo = chr1.compareChromoName(chr2);
				if (compChromo != 0) return compChromo;
			}
		}

		// Start
		if (start > i2.start) return 1;
		if (start < i2.start) return -1;

		// End
		if (end > i2.end) return 1;
		if (end < i2.end) return -1;

		return 0;
	}

	/**
	 * How far apart are these intervals?
	 * @return  Distance or -1 if they are not comparable (i.e. different chromosomes)
	 */
	public int distance(Marker interval) {
		if (!interval.getChromosomeName().equals(getChromosomeName())) return -1;

		if (intersects(interval)) return 0;

		if (start > interval.getEnd()) return start - interval.getEnd();
		if (interval.getStart() > end) return interval.getStart() - end;

		throw new RuntimeException("This should never happen!");
	}

	/**
	 * Go up (parent) until we find an instance of 'clazz'
	 */
	@SuppressWarnings("rawtypes")
	public Interval findParent(Class clazz) {
		if (this.getClass().equals(clazz)) return this;
		if ((parent != null) && (parent instanceof Marker)) return ((Marker) parent).findParent(clazz);
		return null;
	}

	public Chromosome getChromosome() {
		return (Chromosome) findParent(Chromosome.class);
	}

	/**
	 * Find chromosome name
	 * @return
	 */
	public String getChromosomeName() {
		Chromosome chromo = (Chromosome) findParent(Chromosome.class);
		if (chromo != null) return chromo.getId();
		return "";
	}

	/**
	 * Find chromosome and return it's number
	 * 
	 * @return Chromosome number if found, -1 otherwise
	 */
	public double getChromosomeNum() {
		Chromosome chromo = (Chromosome) findParent(Chromosome.class);
		if (chromo != null) return chromo.chromosomeNum;
		return -1;
	}

	/**
	 * Find genome 
	 * @return
	 */
	public Genome getGenome() {
		return (Genome) findParent(Genome.class);
	}

	/**
	 * Find genome name
	 * @return
	 */
	public String getGenomeName() {
		Genome genome = (Genome) findParent(Genome.class);
		if (genome != null) return genome.getId();
		return "";
	}

	@Override
	public Marker getParent() {
		return (Marker) parent;
	}

	public String getType() {
		return type;
	}

	@Override
	public int hashCode() {
		int hashCode = getChromosomeName().hashCode();
		hashCode = hashCode * 31 + start;
		hashCode = hashCode * 31 + end;
		hashCode = hashCode * 31 + strand;
		hashCode = hashCode * 31 + id.hashCode();
		return hashCode;
	}

	/**
	 * Is 'interval' completely included in 'this'?
	 * @return  return true if 'this' includes 'interval' 
	 */
	public boolean includes(Marker interval) {
		if (!interval.getChromosomeName().equals(getChromosomeName())) return false;
		return (start <= interval.start) && (interval.end <= end);
	}

	/**
	 * Do the intervals intersect?
	 * @return  return true if this intersects 'interval' 
	 */
	public boolean intersects(Marker interval) {
		if (!interval.getChromosomeName().equals(getChromosomeName())) return false;
		return (interval.getEnd() >= start) && (interval.getStart() <= end);
	}

	/**
	 * how much do intervals intersect?
	 * @return  number of bases these intervals intersect
	 */
	public int intersectSize(Marker interval) {
		if (!interval.getChromosomeName().equals(getChromosomeName())) return 0;

		int start = Math.max(this.start, interval.getStart());
		int end = Math.min(this.end, interval.getEnd());

		if (end < start) return 0;
		return (end - start) + 1;
	}

	/**
	 * Adjust parent if it does not include child? 
	 * @return
	 */
	protected boolean isAdjustIfParentDoesNotInclude(Marker parent) {
		return false;
	}

	/**
	 * Show an error if parent does not include child? 
	 * @return
	 */
	protected boolean isShowWarningIfParentDoesNotInclude() {
		return false;
	}

	/**
	 * Parse a line (form a file)
	 * Format: "chromosome \t start \t end \t id \n" 
	 */
	public void parse(String line, int lineNum, Genome genome, int positionBase) {
		line = line.trim(); // Remove spaces

		// Ignore empty lines and comment lines
		if ((line.length() > 0) && (!line.startsWith("#"))) {
			// Parse line
			String fields[] = line.split("\\s+");

			// Is line OK?
			if (fields.length >= 3) {
				Chromosome chromo = genome.getChromosome(fields[0].trim());
				if (chromo == null) System.err.println("WARNING: Chromosome '" + fields[0] + "' not found in genome '" + genome.getGenomeName() + "', version '" + genome.getVersion() + "'!\n\tLine: " + lineNum + "\t'" + line + "'");
				parent = chromo;
				start = Gpr.parseIntSafe(fields[1]) - positionBase;
				end = Gpr.parseIntSafe(fields[2]) - positionBase;

				if (fields.length >= 4) {
					// Join all ids using a space character (and remove all spaces)
					for (int t = 3; t < fields.length; t++)
						id += fields[t].trim() + " ";
					id = id.trim();
				}
			} else throw new RuntimeException("Error line " + lineNum + " (number of fields is " + fields.length + "):\t" + line);
		}
	}

	/**
	 * Calculate the effect of this seqChange
	 * @param seqChange
	 * @param changeEffect
	 * @return
	 */
	public List<ChangeEffect> seqChangeEffect(SeqChange seqChange, ChangeEffect changeEffect) {
		if (!intersects(seqChange)) return ChangeEffect.emptyResults(); // Sanity check
		EffectType effType = EffectType.valueOf(type);
		if (effType == null) throw new RuntimeException("Cannot find effect type for '" + type + "'");
		changeEffect.set(this, effType, "");
		return changeEffect.newList();
	}

	@Override
	public String toString() {
		return getChromosomeName() + "\t" + start + "-" + end //
				+ " " //
				+ type + ((id != null) && (id.length() > 0) ? " '" + id + "'" : "");
	}

	@Override
	public String toStringSave() {
		return getChromosomeName() + "\t" + start + "\t" + end + "\t" + id;
	}

}
