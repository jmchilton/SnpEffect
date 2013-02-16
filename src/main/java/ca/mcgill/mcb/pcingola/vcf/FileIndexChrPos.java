package ca.mcgill.mcb.pcingola.vcf;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.channels.FileChannel;
import java.util.HashMap;
import java.util.Set;

import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * Index a file that has "chr \t pos" as the beginning of a line (e.g. VCF)
 * 
 * WARNING: It is assumed that the file is ordered by position (chromosome order does not matter)
 * 
 * TODO: Rewrite this method using java.io.RandomAccessFile
 * 
 * @author pcingola
 */
public class FileIndexChrPos {

	/**
	 * A part of a file
	 * @author pcingola
	 *
	 */
	class FileRegion {
		long start, end;
		String lineStart, lineEnd;
	}

	/**
	 * A line and the position on the file where it begins
	 * @author pcingola
	 *
	 */
	class LineAndPos {
		String line;
		long position;

		@Override
		public String toString() {
			String str = "";
			if (line != null) {
				if (line.length() > 50) str = line.substring(0, 49) + "...";
				else str = line;
			}
			return position + "\t" + str;
		}
	}

	public static final int POS_OFFSET = 1; // VCF files are one-based
	private static final int BUFF_SIZE = 1024 * 1024;

	boolean verbose = false;
	boolean debug = false;
	String fileName;
	long size = 0;
	RandomAccessFile file;
	HashMap<String, FileRegion> fileRegions = new HashMap<String, FileIndexChrPos.FileRegion>(); // Store file regions by chromosome

	public FileIndexChrPos(String fileName) {
		this.fileName = fileName;
	}

	/**
	 * Get chromosome info
	 * @param line
	 * @return
	 */
	String chromo(String line) {
		if (line.startsWith("#")) return null;
		return line.split("\\t")[0];
	}

	/**
	 * Close file
	 */
	public void close() {
		try {
			if (file != null) file.close();
		} catch (IOException e) {
			System.err.println("I/O problem while closing file '" + fileName + "'");
			throw new RuntimeException(e);
		}
		file = null;
	}

	/**
	 * Dump a region of the file to STDOUT
	 * @param start
	 * @param end
	 */
	void dump(long start, long end) {
		if (verbose) System.err.println("\tDumping file '" + fileName + "' interval [ " + start + " , " + end + " ]");
		try {
			byte buff[] = new byte[BUFF_SIZE];
			file.seek(start);
			for (long curr = start; curr <= end;) {
				long len = Math.min(BUFF_SIZE, end - curr + 1); // Maximum length to read
				int read = file.read(buff, 0, (int) len); // Read file

				if (read <= 0) break; // Error or nothing read, abort

				String out = new String(buff, 0, read);
				System.out.print(out);

				curr += read;
			}
		} catch (Exception e) {
			throw new RuntimeException("Error reading file '" + fileName + "' from position " + start + " to " + end);
		}
	}

	/**
	 * Dump all lines in the interval chr:posStart-posEnd
	 * 
	 * @param chr
	 * @param posStart
	 * @param posEnd
	 */
	public void dump(String chr, int posStart, int posEnd) {
		debug = true;
		Gpr.debug("DEBUG!");
		long fileStart = find(chr, posStart, false);
		long fileEnd = find(chr, posEnd, true);

		dump(fileStart, fileEnd - 1);
	}

	/**
	 * Find the position in the file for the first character of the first line whose genomic position is less or equal than 'chrPos'
	 * @param chrPos
	 * @param start
	 * @param lineStart
	 * @param end
	 * @param lineEnd
	 * @return
	 */
	long find(int chrPos, long start, String lineStart, long end, String lineEnd, boolean endOfLine) {
		//---
		// Check break conditions
		//---

		// Is it lineStart?
		int posStart = pos(lineStart);
		if (chrPos <= posStart) {
			if (endOfLine) return start + lineStart.length(); // End of line
			else return start; // Begining of 'start' line
		}

		// Is it lineEnd?
		int posEnd = pos(lineEnd);
		if (posEnd <= chrPos) {
			if (endOfLine) return end + lineEnd.length(); // End of line
			else return end; // Begining of 'start' line
		}

		// Interval is less than one line?
		if ((start + lineStart.length()) >= end) {
			if (endOfLine) return start + lineStart.length(); // End of line
			else return start; // Begining of 'start' line
		}

		if (debug) Gpr.debug("Find:\t" + chrPos + "\t[" + posStart + ", " + posEnd + "]\tFile: [" + start + " , " + end + "]\tsize: " + (end - start));

		// Sanity check
		if (posStart >= posEnd) throw new RuntimeException("This should never happen! Is the file sorted by position?");

		//---
		// Recurse
		//---
		long mid = (start + end) / 2;
		String lineMid = getLine(mid).line;
		long posMid = pos(lineMid);

		if (chrPos <= posMid) return find(chrPos, start, lineStart, mid, lineMid, endOfLine);
		else return find(chrPos, mid, lineMid, end, lineEnd, endOfLine);
	}

	/**
	 * Find the position in the file for the first character of the first line equal or less than a specific chr:pos
	 * @param chr
	 * @param pos
	 * @return
	 */
	long find(String chr, int pos, boolean lessEq) {
		chr = Chromosome.simpleName(chr);
		FileRegion fr = fileRegions.get(chr);
		if (fr == null) throw new RuntimeException("No such chromosome: '" + chr + "'");
		long posFound = find(pos, fr.start, fr.lineStart, fr.end, fr.lineEnd, lessEq);
		return getLine(posFound).position;
	}

	/**
	 * Get a byte from a file
	 * @param bytePosition
	 * @return
	 */
	public byte get(long bytePosition) {
		try {
			// Change position if needed
			if (file.getFilePointer() != bytePosition) file.seek(bytePosition);
			return (byte) file.read();
		} catch (IOException e) {
			throw new RuntimeException("Error readin file '" + fileName + "' at position " + bytePosition, e);
		}
	}

	public byte[] get(long bytePosition, int len) {
		try {
			byte buff[] = new byte[len];

			// Change position if needed
			if (file.getFilePointer() != bytePosition) file.seek(bytePosition);

			int read = file.read(buff);

			// Nothing to read?
			if (read <= 0) return null;

			// Buffer was too long? Return an array of byte with exactly the number of byte that were  
			if (read < buff.length) {
				byte newBuff[] = new byte[read];
				System.arraycopy(buff, 0, newBuff, 0, read);
				buff = newBuff;
			}

			return buff;
		} catch (IOException e) {
			throw new RuntimeException("Error readin file '" + fileName + "' at position " + bytePosition, e);
		}
	}

	/**
	 * Available chromosomes
	 * @return
	 */
	public Set<String> getChromos() {
		return fileRegions.keySet();
	}

	/**
	 * Get position where 'chr' ends
	 * @param chr
	 * @return -1 if 'chr' is not in the index
	 */
	public long getEnd(String chr) {
		chr = Chromosome.simpleName(chr);
		FileRegion fr = fileRegions.get(chr);
		if (fr == null) return -1;
		return fr.end;
	}

	/**
	 * Get file region for a given chrosmome
	 * @param chr
	 * @return
	 */
	FileRegion getFileRegion(String chr) {
		chr = Chromosome.simpleName(chr);
		FileRegion fr = fileRegions.get(chr);
		if (fr == null) {
			fr = new FileRegion();
			fileRegions.put(chr, fr);
		}
		return fr;
	}

	/**
	 * Get the line where 'pos' hits
	 * 
	 * TODO: This is really slow for huge files and huge lines. I should optimize this.
	 * 
	 * @param pos
	 * @return A string with the line that 'pos' hits, null if it's out of boundaries
	 */
	public LineAndPos getLine(long pos) {
		long size = size();
		if ((pos >= size) || (pos < 0)) return null;

		LineAndPos linePos = new LineAndPos();
		StringBuffer sb = new StringBuffer();

		// Get bytes before 'pos'
		long position;
		for (position = pos - 1; position >= 0; position--) {
			byte b = get(position);
			if (b == '\n') break;
			sb.insert(0, (char) b);
		}
		linePos.position = position + 1;

		// Get bytes after 'pos'
		for (position = pos; position < size; position++) {
			byte b = get(position);
			if (b == '\n') break;
			sb.append((char) b);
		}
		linePos.line = sb.toString();

		if (debug) Gpr.debug("Line & Position: " + linePos);
		return linePos;
	}

	/**
	 * Get position where 'chr' starts
	 * @param chr
	 * @return -1 if 'chr' is not in the index
	 */
	public long getStart(String chr) {
		chr = Chromosome.simpleName(chr);
		FileRegion fr = fileRegions.get(chr);
		if (fr == null) return -1;
		return fr.start;
	}

	/**
	 * Index chromosomes in the whole file 
	 */
	public void index() {
		// Last line (minus '\n' character, minus one)
		long end = size() - 1;
		String lineEnd = getLine(end).line;
		String chrEnd = chromo(lineEnd);

		// Add fileRegion.end for last chromsome in the file
		FileRegion fr = getFileRegion(chrEnd);
		fr.end = end;
		fr.lineEnd = lineEnd;
		if (verbose) System.err.println("\tindex:\t" + chrEnd + "\t" + end);

		// Find first non-comment line
		long start = 0;
		String lineStart = "";
		for (start = 0; start < size; start += lineStart.length() + 1) {
			lineStart = getLine(start).line;
			if (chromo(lineStart) != null) break;
		}

		String chrStart = chromo(lineStart);

		// Add fileRegion.start for first chromsome in the file
		fr = getFileRegion(chrStart);
		fr.start = start;
		fr.lineStart = lineStart;
		if (verbose) System.err.println("\tindex:\t" + chrStart + "\t" + start);

		// Index the rest of the file
		indexChromos(start, lineStart, end, lineEnd);
	}

	/**
	 * Index chromosomes in a region of a file
	 * @param start
	 * @param lineStart
	 * @param end
	 * @param lineEnd
	 */
	void indexChromos(long start, String lineStart, long end, String lineEnd) {
		if (debug) Gpr.debug("Index:"//
				+ "\n\t" + start + "(" + (((double) start) / size()) + ") :\t" + s(lineStart) //
				+ "\n\t" + end + "(" + (((double) end) / size()) + ") :\t" + s(lineEnd));

		if (start > end) throw new RuntimeException("This should never happen! Start: " + start + "\tEnd: " + end);

		String chrStart = chromo(lineStart);
		String chrEnd = chromo(lineEnd);

		if (chrStart.equals(chrEnd)) {
			if (debug) Gpr.debug("Chromo:\tlineStart: " + chrStart + "\tlineEnd: " + chrEnd + "\t==> Back!");
			return;
		}
		if (debug) Gpr.debug("Chromo:\tlineStart: " + chrStart + "\tlineEnd: " + chrEnd);

		if ((start + lineStart.length() + 1) >= end) {
			if (verbose) System.err.println("\t\t" + chrStart + " / " + chrEnd + "\t" + start + " / " + end);

			// Add index where chromosome starts
			getFileRegion(chrEnd).start = getLine(end).position;
			getFileRegion(chrEnd).lineStart = lineEnd;

			// Add index where chromosome ends
			getFileRegion(chrStart).end = getLine(start).position;
			getFileRegion(chrStart).lineEnd = lineStart;
			return;
		}

		long mid = (start + end) / 2;
		String lineMid = getLine(mid).line;
		if (debug) Gpr.debug("Mid: " + mid + "\t" + s(lineMid));

		if (debug) Gpr.debug("First half recustion:");
		indexChromos(start, lineStart, mid, lineMid);

		if (debug) Gpr.debug("Second half recustion:");
		indexChromos(mid, lineMid, end, lineEnd);
	}

	void init(FileChannel channel) throws IOException {
	}

	/**
	 * Open file and initiate mappings
	 */
	public void open() {
		try {
			File f = new File(fileName);
			size = f.length();
			file = new RandomAccessFile(f, "r");
		} catch (FileNotFoundException e) {
			System.err.println("File not found '" + fileName + "'");
			throw new RuntimeException(e);
		}
	}

	/**
	 * The position argument of a line, or zero if not found
	 * @param line
	 * @return
	 */
	int pos(String line) {
		if (line.startsWith("#")) return 0; // In VCF, positions are one-based, so zero denotes an error
		return Gpr.parseIntSafe(line.split("\\t")[1]) - POS_OFFSET;
	}

	String s(String s) {
		if (s == null) return "null";
		return s.length() <= 50 ? s : s.substring(0, 50) + "...";
	}

	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}

	/**
	 * File size
	 * @return
	 */
	public long size() {
		return size;
	}
}
