package ca.mcgill.mcb.pcingola.interval;

import java.io.Serializable;

import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * A genomic interval.
 * Note: Intervals are assumed to be zero-based and inclusive
 *       i.e. an interval including the first base up to base X would 
 *       be [0,X] NOT [1,X] 
 * 
 * @author pcingola
 */
public class Interval implements Comparable<Interval>, Serializable {

	private static final long serialVersionUID = -3547434510230920403L;

	protected int start, end;
	protected int strand;
	protected String id = ""; // Interval's ID (e.g. gene name, transcript ID)
	protected Interval parent;

	public Interval(Interval parent, int start, int end, int strand, String id) {
		if (start > end) throw new RuntimeException("Interval error: end before start. Start:" + start + ", End: " + end);
		this.start = start;
		this.end = end;
		this.id = id;
		this.strand = strand;
		this.parent = parent;
	}

	/** 
	 * Parse a line (form a file)
	 * Format: "chromosome \t start \t end \t id \n"
	 */
	public Interval(String line, int lineNum) {
		parse(line, lineNum);
	}

	/**
	 * Compare by start and end
	 */
	@Override
	public int compareTo(Interval i2) {
		// Start
		if (start > i2.start) return 1;
		if (start < i2.start) return -1;

		// End
		if (end > i2.end) return 1;
		if (end < i2.end) return -1;

		return 0;
	}

	public boolean equals(Interval interval) {
		return compareTo(interval) == 0;
	}

	public int getEnd() {
		return end;
	}

	public String getId() {
		return id;
	}

	public Interval getParent() {
		return parent;
	}

	public int getStart() {
		return start;
	}

	public int getStrand() {
		return strand;
	}

	@Override
	public int hashCode() {
		int hashCode = 0;
		hashCode = hashCode * 31 + start;
		hashCode = hashCode * 31 + end;
		hashCode = hashCode * 31 + strand;
		hashCode = hashCode * 31 + id.hashCode();
		return hashCode;
	}

	/**
	 * @return  return true if this intersects 'interval' 
	 */
	public boolean intersects(Interval interval) {
		return (interval.getEnd() >= start) && (interval.getStart() <= end);
	}

	/**
	 * @return  true if this interval contains point (inclusive)
	 */
	public boolean intersects(long point) {
		return (start <= point) && (point <= end);
	}

	public boolean isStrandMinus() {
		return strand < 0;
	}

	public boolean isStrandPlus() {
		return strand >= 0;
	}

	public boolean isValid() {
		return (start >= 0) && (start <= end);
	}

	/**
	 * Parse a line (form a file)
	 * Format: "start \t end \t id \n" 
	 * 
	 * @param line
	 * @param lineNum
	 */
	void parse(String line, int lineNum) {
		line = line.trim(); // Remove spaces

		// Ignore empty lines and comment lines
		if ((line.length() > 0) && (!line.startsWith("#"))) {
			// Parse line
			String fields[] = line.split("\\s+");

			// Is line OK?
			if (fields.length >= 2) {
				start = Gpr.parseIntSafe(fields[0]);
				end = Gpr.parseIntSafe(fields[1]);

				if (fields.length >= 3) {
					// Join all ids using a space character (and remove all spaces)
					for (int t = 2; t < fields.length; t++)
						id += fields[t].trim() + " ";
					id = id.trim();
				}
			} else throw new RuntimeException("Error line " + lineNum + " (number of fields is " + fields.length + "):\t" + line);
		}
	}

	public void setEnd(int end) {
		this.end = end;
	}

	public void setId(String id) {
		this.id = id;
	}

	public void setStart(int start) {
		this.start = start;
	}

	public void setStrand(int strand) {
		this.strand = strand;
	}

	public int size() {
		return end - start + 1;
	}

	@Override
	public String toString() {
		return start + "-" + end //
				+ ((id != null) && (id.length() > 0) ? " '" + id + "'" : "");
	}

	/**
	 * Show it as an ASCII art
	 * @param maxLen
	 * @return
	 */
	public String toStringAsciiArt(int maxLen) {
		StringBuilder sb = new StringBuilder();

		for (int i = 0; i < maxLen; i++) {
			if ((i >= start) && (i <= end)) sb.append('-');
			else sb.append(' ');
		}

		return sb.toString();
	}

	public String toStringSave() {
		return start + "\t" + end + "\t" + id;
	}

}
