package ca.mcgill.mcb.pcingola.snpEffect.factory;

import ca.mcgill.mcb.pcingola.genBank.Feature;
import ca.mcgill.mcb.pcingola.genBank.Feature.Type;
import ca.mcgill.mcb.pcingola.genBank.Features;
import ca.mcgill.mcb.pcingola.interval.Cds;
import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.SnpEffectPredictor;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.GprSeq;

/**
 * This class creates a SnpEffectPredictor from a 'features' file.
 * This includes derived formats as GenBank and Embl.
 * 
 * References:
 * 		http://www.ebi.ac.uk/embl/Documentation/User_manual/printable.html
 * 		http://www.ebi.ac.uk/embl/Documentation/FT_definitions/feature_table.html
 * 
 * 
 * @author pcingola
 */
public abstract class SnpEffPredictorFactoryFeatures extends SnpEffPredictorFactory {

	public static final int OFFSET = 1;
	Chromosome chromosome; // It is assumed that there is only one 'Chromosome' (i.e. only one 'SOURCE' feature)

	public SnpEffPredictorFactoryFeatures(Config config) {
		super(config, OFFSET);
	}

	/**
	 *	Add all features 
	 */
	protected void addFeatures(Features features) {
		//---
		// Add all chromosome 
		//---
		for (Feature f : features.getFeatures()) {
			// Convert coordinates to zero-based 
			int start = f.getStart() - inOffset;
			int end = f.getEnd() - inOffset;

			// Add chromosome
			if (f.getType() == Type.SOURCE) {
				if (chromosome != null) throw new RuntimeException("SOURCE already assigned to chromosome");
				chromosome = new Chromosome(genome, start, end, 1, chromoName(f)); // Convert coordinates to zero-based 
				add(chromosome);
			}
		}

		// Sanity check
		if (chromosome == null) throw new RuntimeException("Could not find SOURCE feature");

		//---
		// Add a genes
		//---
		for (Feature f : features.getFeatures()) {
			if (f.getType() == Type.GENE) findOrCreateGene(f, chromosome, false);
		}

		//---
		// Add transcripts
		//---
		for (Feature f : features.getFeatures()) {
			// Transcript
			if (f.getType() == Type.MRNA) {
				// Convert coordinates to zero-based 
				int start = f.getStart() - inOffset;
				int end = f.getEnd() - inOffset;

				String trId = getTrId(f);
				Gene gene = findOrCreateGene(f, chromosome, false); // Find or create gene

				// Add transcript
				Transcript tr = new Transcript(gene, start, end, f.isComplement() ? 1 : -1, trId);
				add(tr);
			}
		}

		//---
		// Add CDS and protein coding information 
		//---
		for (Feature f : features.getFeatures()) {
			// Add CDS
			if (f.getType() == Type.CDS) {
				// Convert coordinates to zero-based 
				int start = f.getStart() - inOffset;
				int end = f.getEnd() - inOffset;

				// Try to fing transcript
				String trId = getTrId(f);

				Transcript tr = findTranscript(trId);
				if (tr == null) {
					if (verbose) System.err.println("WARNING: Transcript '" + trId + "' not found. Creating new transcript." + f);

					// Not found? => Create gene and transcript
					Gene gene = findOrCreateGene(f, chromosome, false); // Find or create gene
					trId = "Tr_" + start + "_" + end;
					tr = new Transcript(gene, start, end, f.isComplement() ? 1 : -1, trId);
					add(tr);
				}

				// Mark transcript as protein coding
				if (f.get("translation") != null) tr.setProteinCoding(true);

				Cds cds = new Cds(tr, f.getStart() - inOffset, f.getEnd() - inOffset, f.isComplement() ? 1 : -1, "CDS_" + trId);
				add(cds);
			}
		}
	}

	/**
	 * Find or create a chromosome name for a feature
	 * @param f
	 * @return
	 */
	String chromoName(Feature f) {
		if (f.getType() != Type.SOURCE) throw new RuntimeException("Cannot find chromosome name in a non-SOURCE feature");
		String chrName = f.get("chromosome");
		return chrName != null ? chrName : "chr1";
	}

	@Override
	public SnpEffectPredictor create() {
		// Read gene intervals from a file
		System.out.println("Reading data file  : '" + fileName + "'");
		try {
			// Read file and add all features
			Features features = readFeatures();
			addFeatures(features);

			// Some clean-up before readng exon sequences
			beforeExonSequences();

			// Get exon sequences
			String sequence = sequence(features);
			addExonSequences(chromosome.getId(), sequence);

			// Finish up (fix problems, add missing info, etc.)
			finishUp(false);

			// Check that exons have sequences
			boolean error = config.getGenome().showStats();
			if (error) throw new RuntimeException("Most Exons do not have sequences!");
		} catch (Exception e) {
			e.printStackTrace();
			throw new RuntimeException("Error reading file '" + fileName + "'\n" + e);
		}

		return snpEffectPredictor;
	}

	//	/**
	//	 * Set exon sequences
	//	 * @param sequence
	//	 */
	//	protected void exonSequences(String sequence) {
	//		for (Gene g : genome.getGenes())
	//			for (Transcript t : g)
	//				for (Exon e : t) {
	//					if (e.getEnd() < sequence.length()) {
	//						String seq = sequence.substring(e.getStart(), e.getEnd() + 1);
	//						e.setSequence(seq);
	//					} else throw new RuntimeException("WARNING: Exon outside sequence. Sequence length is " + sequence.length() + ", exon range is [ " + e.getStart() + " , " + e.getEnd() + " ]");
	//				}
	//	}

	Gene findOrCreateGene(Feature f, Chromosome chr, boolean warn) {
		int start = f.getStart() - inOffset;
		int end = f.getEnd() - inOffset;

		String geneId = geneId(f, start, end);

		Gene gene = findGene(geneId);
		if (gene == null) {
			gene = new Gene(chr, start, end, f.isComplement() ? 1 : -1, geneId, geneId, "");
			add(gene);
			if (warn) System.err.println("WARNING: Gene '" + geneId + "' not found");
		}

		return gene;
	}

	/**
	 * Try to get geneIDs
	 * @param f
	 * @return
	 */
	protected String geneId(Feature f, int start, int end) {
		// Try 'gene'...
		String geneId = f.get("gene");
		if (geneId != null) return geneId;

		// Try 'locus'...
		geneId = f.get("locus_tag");
		if (geneId != null) return geneId;

		// Try 'locus'...
		geneId = f.get("locus_tag");
		if (geneId != null) return geneId;

		return "Gene_" + start + "_" + end;
	}

	/**
	 * Create a transciript ID based on a feature
	 * @param f
	 * @return
	 */
	protected String getTrId(Feature f) {
		String id = f.get("product");
		if (id != null) id = id.replaceAll("\\s", "_");
		if ((id == null) || (id.length() > 20) && (f.get("gene") != null)) return "Tr_" + f.get("gene");
		return id;
	}

	/**
	 * Get features from file
	 * @return
	 */
	protected abstract Features readFeatures();

	/**
	 * Get sequence either from features or from FASTA file
	 * @param features
	 * @return
	 */
	String sequence(Features features) {
		String seq = features.getSequence();
		if ((seq != null) && !seq.isEmpty()) return seq;
		if (verbose) System.out.println("No sequence found in feature file.");

		// No sequence information in 'features' file? => Try to read a sequence from a fasta file
		for (String fastaFile : config.getFileNameGenomeFasta()) {
			if (verbose) System.out.println("\tTrying fasta file '" + fastaFile + "'");

			if (Gpr.canRead(fastaFile)) {
				seq = GprSeq.fastaSimpleRead(fastaFile);
				if ((seq != null) && !seq.isEmpty()) return seq;
			}
		}

		throw new RuntimeException("Cannot find sequence for '" + config.getGenome().getVersion() + "'");
	}
}
