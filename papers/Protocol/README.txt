
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



