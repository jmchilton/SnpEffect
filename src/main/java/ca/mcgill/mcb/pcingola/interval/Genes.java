package ca.mcgill.mcb.pcingola.interval;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

/**
 * A collection of genes (marker intervals)
 * Note: It is assumed that all genes belong to the same genome
 * 
 * @author pcingola
 */
@SuppressWarnings("serial")
public class Genes implements Iterable<Gene>, Serializable {

	public boolean debug = false;
	Genome genome;
	HashMap<String, Gene> genesById;

	public Genes(Genome genome) {
		genesById = new HashMap<String, Gene>();
		this.genome = genome;
	}

	/**
	 * Add a gene interval to this collection
	 * @param gene
	 */
	public void add(Gene gene) {
		genesById.put(gene.getId(), gene);
	}

	/**
	 * Creates a list of splice sites intervals
	 * For a definition of splice site, see comments at the beginning of SpliceSite.java
	 */
	public Collection<Marker> createSpliceSites() {
		HashMap<String, Marker> map = new HashMap<String, Marker>(); // Use a map in order to remove repeated splice sites (different transcripts may have the same exons)

		// For each gene, transcript and exon
		for( Gene gene : this ) {
			for( Transcript tr : gene ) {
				List<SpliceSite> slist = tr.createSpliceSites();

				// Store all markers in hash
				for( SpliceSite ss : slist ) {
					String key = ss.getClass().getSimpleName() + " " + ss.getChromosomeName() + ":" + ss.getStart() + "-" + ss.getEnd() + "_" + ss.getId();
					map.put(key, ss);
				}
			}
		}

		return map.values();
	}

	/**
	 * Creates a list of UP/DOWN stream regions (for each transcript)
	 * Upstream (downstream) stream is defined as upDownLength before (after) transcript
	 * 
	 * Note: If upDownLength <=0 no interval is created
	 */
	public List<Marker> createUpDownStream(int upDownLength) {
		ArrayList<Marker> list = new ArrayList<Marker>();
		if( upDownLength <= 0 ) return list;

		// For each gene, transcript
		for( Gene gene : this ) {
			for( Transcript tr : gene ) {
				tr.createUpDownStream(upDownLength);
				list.add(tr.getUpstream());
				list.add(tr.getDownstream());
			}
		}
		return list;
	}

	/**
	 * Obtain a gene interval
	 * @param geneId
	 * @return
	 */
	public Gene get(String geneId) {
		return genesById.get(geneId);
	}

	@Override
	public Iterator<Gene> iterator() {
		return genesById.values().iterator();
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for( Gene gint : this )
			sb.append(gint + "\n");
		return sb.toString();
	}

	public Collection<Gene> values() {
		return genesById.values();
	}
}
