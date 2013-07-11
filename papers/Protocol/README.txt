
Create an "artificial" mutation:
	- Homozygous in affected childs
	- Het in 2 parents
	- Het in two grand-parents

Selected mutation: 
	- Gene CFTR (cystic fibrosis transmembrane conductance regulator)
	- G542*

	- ClinVar: http://www.ncbi.nlm.nih.gov/clinvar/RCV000007535/
		$ grep RCV000007535 clinvar_00-latest.vcf
		7	117227832	rs113993959	G	T	.	.	RS=113993959;RSPOS=117227832;dbSNPBuildID=132;SSR=0;SAO=1;VP=0x050268000601040002110100;GENEINFO=CFTR:1080;WGT=1;VC=SNV;PM;PMC;S3D;NSN;REF;VLD;OTHERKG;LSD;OM;CLNALLE=1;CLNHGVS=NC_000007.13:g.117227832G>T;CLNSRC=GTR|OMIM Allelic Variant|OMIM Allelic Variant;CLNORIGIN=1;CLNSRCID=GTR000500233|602421.0009|602421.0095;CLNSIG=5;CLNDSDB=GeneReviews:NCBI:OMIM:Orphanet:SNOMED CT;CLNDSDBID=NBK1250:C0010674:219700:586:190905008;CLNDBN=Cystic fibrosis;CLNACC=RCV000007535.1

		$ grep RCV000007535 clinvar_00-latest.vcf | cut -f 8 | tr ";" "\n"
		RS=113993959
		RSPOS=117227832
		dbSNPBuildID=132
		SSR=0
		SAO=1
		VP=0x050268000601040002110100
		GENEINFO=CFTR:1080
		WGT=1
		VC=SNV
		PM
		PMC
		S3D
		NSN
		REF
		VLD
		OTHERKG
		LSD
		OM
		CLNALLE=1
		CLNHGVS=NC_000007.13:g.117227832G>T
		CLNSRC=GTR|OMIM Allelic Variant|OMIM Allelic Variant
		CLNORIGIN=1
		CLNSRCID=GTR000500233|602421.0009|602421.0095
		CLNSIG=5
		CLNDSDB=GeneReviews:NCBI:OMIM:Orphanet:SNOMED CT
		CLNDSDBID=NBK1250:C0010674:219700:586:190905008
		CLNDBN=Cystic fibrosis
		CLNACC=RCV000007535.1

Added in samples:
		$ grep -v "^#" mutation.vcf | cut -f 10-
		1/0	1/0	1/1	.	.	1/0	1/0	1/0	1/1	1/1	.	1/0	1/0	.	.	1/0	.

----------------------------------------------------------------------------------------------------



----------------------------------------------------------------------------------------------------

Protocol
--------

Step 1: Download SnpEff , Unzip

Step 2: Download GRCh37.71

Step 3: Download dataset

Step 4: Annotate 

		snpeff -v -noStats -lof -motif -nextProt $genome $vcf > $vcfout

		Notice '-lof' command line option, 
		- We could ommit -nextProt

Step 5: Filter homozygous on cases

		snpsift caseControl -v -tfam Utah-Pedigree-1463.tfam by_chromo_CEPH/eff/cg_panel.vcf > cg_panel.cc.eff.vcf

Step 6: Filter high impact coding variants

		cat cg_panel.cc.eff.vcf | snpsift filter "(Cases[0] = 3) & (Controls[0] = 0) & (EFF[*].IMPACT = 'HIGH')" > filtered.vcf

Step 7: Show effects

		$ grep -v "^#" filtered.vcf | cut -f 8 | tr ";" "\n" | grep "^EFF=" | tr "," "\n"
		EFF=DOWNSTREAM(MODIFIER|||||CFTR|processed_transcript|CODING|ENST00000472848||1)
		NEXT_PROT[beta_strand](LOW||||1419|CFTR|protein_coding|CODING|ENST00000454343|11|1)
		NEXT_PROT[beta_strand](LOW||||1480|CFTR|protein_coding|CODING|ENST00000003084|12|1)
		NEXT_PROT[domain:ABC_transporter_1](LOW||||1480|CFTR|protein_coding|CODING|ENST00000003084|10|1)
		NEXT_PROT[domain:ABC_transporter_1](LOW||||1480|CFTR|protein_coding|CODING|ENST00000003084|11|1)
		NEXT_PROT[domain:ABC_transporter_1](LOW||||1480|CFTR|protein_coding|CODING|ENST00000003084|12|1)
		NEXT_PROT[domain:ABC_transporter_1](LOW||||1480|CFTR|protein_coding|CODING|ENST00000003084|13|1)
		NEXT_PROT[domain:ABC_transporter_1](LOW||||1480|CFTR|protein_coding|CODING|ENST00000003084|14|1)
		NEXT_PROT[topological_domain:Cytoplasmic](LOW||||1480|CFTR|protein_coding|CODING|ENST00000003084|10|1)
		NEXT_PROT[topological_domain:Cytoplasmic](LOW||||1480|CFTR|protein_coding|CODING|ENST00000003084|11|1)
		NEXT_PROT[topological_domain:Cytoplasmic](LOW||||1480|CFTR|protein_coding|CODING|ENST00000003084|12|1)
		NEXT_PROT[topological_domain:Cytoplasmic](LOW||||1480|CFTR|protein_coding|CODING|ENST00000003084|13|1)
		NEXT_PROT[topological_domain:Cytoplasmic](LOW||||1480|CFTR|protein_coding|CODING|ENST00000003084|14|1)
		NEXT_PROT[topological_domain:Cytoplasmic](LOW||||1480|CFTR|protein_coding|CODING|ENST00000003084|15|1)
		NEXT_PROT[topological_domain:Cytoplasmic](LOW||||1480|CFTR|protein_coding|CODING|ENST00000003084|8|1)
		NEXT_PROT[topological_domain:Cytoplasmic](LOW||||1480|CFTR|protein_coding|CODING|ENST00000003084|9|1)
		STOP_GAINED(HIGH|NONSENSE|Gga/Tga|G481*|1419|CFTR|protein_coding|CODING|ENST00000454343|11|1)
		STOP_GAINED(HIGH|NONSENSE|Gga/Tga|G512*|1437|CFTR|protein_coding|CODING|ENST00000426809|11|1|WARNING_TRANSCRIPT_INCOMPLETE)
		STOP_GAINED(HIGH|NONSENSE|Gga/Tga|G542*|1480|CFTR|protein_coding|CODING|ENST00000003084|12|1)
		UPSTREAM(MODIFIER|||||AC000111.5|processed_pseudogene|NON_CODING|ENST00000448200||1)

Step 8 : Show pedigree		
		$ mkdir pedShow
		[eq8302@ehs cg_panel]$ snpsift pedShow Utah-Pedigree-1463.tfam filtered.vcf pedShow/
		00:00:00.000	Reading vcf file 'filtered.vcf'
		00:00:00.049	Drawing pedigree for '7:117227832', output dir: pedShow//7_117227832
		00:00:00.058	Done

