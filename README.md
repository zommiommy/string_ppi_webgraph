# string_ppi_webgraph
Experiments on merging and compressing the whole inter-species string ppi graph using webgraph-rs

# Sources
Needed files from [string ppi](https://string-db.org/cgi/download?sessionId=btvSejFLpDnW):
* [`protein.links.full.v12.0.txt.gz`](https://stringdb-downloads.org/download/protein.links.full.v12.0.txt.gz) 200GB
```
protein1 protein2 neighborhood neighborhood_transferred fusion cooccurence homology coexpression coexpression_transferred experiments experiments_transferred database database_transferred textmining textmining_transferred combined_score
23.BEL05_00025 23.BEL05_01210 0 0 0 153 0 0 0 0 0 0 0 0 0 153
23.BEL05_00025 23.BEL05_13800 0 0 0 167 0 0 0 0 0 0 0 0 0 167
23.BEL05_00025 23.BEL05_13205 0 91 0 0 0 0 45 0 0 0 126 0 61 192
23.BEL05_00025 23.BEL05_12585 0 56 0 0 0 0 267 0 0 0 0 0 66 297
23.BEL05_00025 23.BEL05_17525 0 0 0 181 0 0 167 0 0 0 0 0 0 288
23.BEL05_00025 23.BEL05_12855 0 0 0 0 0 0 57 0 53 0 181 0 57 218
23.BEL05_00025 23.BEL05_18145 0 0 0 179 0 0 0 0 0 0 0 0 0 179
23.BEL05_00025 23.BEL05_04580 0 0 0 180 0 0 0 0 0 0 0 0 0 180
23.BEL05_00025 23.BEL05_10565 0 55 0 0 0 0 109 0 0 0 184 0 56 264
```
* [`enrichment.terms.v12.0`](https://stringdb-downloads.org/download/protein.enrichment.terms.v12.0.txt.gz) 24GB
```
#string_protein_id	category	term	description
100053.GCA_002009845_00001	Biological Process (Gene Ontology)	GO:0006139	Nucleobase-containing compound metabolic process
100053.GCA_002009845_00001	Biological Process (Gene Ontology)	GO:0006396	RNA processing
100053.GCA_002009845_00001	Biological Process (Gene Ontology)	GO:0006399	tRNA metabolic process
100053.GCA_002009845_00001	Biological Process (Gene Ontology)	GO:0006725	Cellular aromatic compound metabolic process
100053.GCA_002009845_00001	Biological Process (Gene Ontology)	GO:0006807	Nitrogen compound metabolic process
100053.GCA_002009845_00001	Biological Process (Gene Ontology)	GO:0008033	tRNA processing
100053.GCA_002009845_00001	Biological Process (Gene Ontology)	GO:0008152	Metabolic process
100053.GCA_002009845_00001	Biological Process (Gene Ontology)	GO:0008616	Queuosine biosynthetic process
100053.GCA_002009845_00001	Biological Process (Gene Ontology)	GO:0009058	Biosynthetic process
```
* [`protein.aliases.v12.0.txt.gz`]()
```
#string_protein_id	alias	source
23.BEL05_00025	23	RefSeq_xref_taxon
23.BEL05_00025	A0A1E5IUA8	UniProt_AC
23.BEL05_00025	A0A1E5IUA8_SHECO	UniProt_ID
23.BEL05_00025	BEL05_00025	RefSeq
23.BEL05_00025	BEL05_00025	RefSeq_locus
23.BEL05_00025	BEL05_00025	UniProt_GN_ORFNames
23.BEL05_00025	FAD-dependent oxidoreductase	UniProt_DE_SubName_Full
23.BEL05_00025	Hypothetical protein	RefSeq_product
23.BEL05_00025	NZ_MCBT01000037.1	UniProt_DR_RefSeq
23.BEL05_00025	WP_069671263.1	UniProt_DR_RefSeq
23.BEL05_00025	hypothetical protein	RefSeq_product
```
* [`oma-groups`](https://omabrowser.org/All/oma-uniprot.txt.gz)  has the cliques with OMA IDs
```
# Orthologous groups from OMA release of All.Jun2023
# This release has 1251567 groups covering 15258796 proteins from 2851 species
# Format: group number<tab>Fingerprint<tab>tab-separated list of OMA Entry IDs
1	FDRGWTQ	HALJB02176	HALHT01804	HALMA01103	HALMD02105	HALUD02830	NATM801193	NATPD01971	HALMT03210	HALVU02175
2	NDVCYKD	UC00500941	UC01800155
3	ERMATNI	HALHT02314	HALMD03006	HALUD00015	NATM800876	NATPD02565	HALMT00683	HALVD00730	HALBP02408	HALWD00195	HALWC00213	HALLT00193	HALVU02374	HALXS02236	NATA101400	NATP100779	NATGS01263
4	IGEKYQC	ARCFU00272	FERPA02377
5	PEPLMYV	PYRAE01670	PYRAR01391	PYRCJ00325	PYRIL00858	PYRNV01190	PYROT00738	THETK00933	THEU700615
6	QCVFVDN	UC00500059	ARCG500944
7	LECAECK	PYRAE01463	PYRAR02059	PYRCJ01039	PYRIL00946	PYRNV01296	PYROT02751	THETK01862	THEU701530
```
* [`oma-species`](https://omabrowser.org/All/oma-species.txt)
```
# Mapping of OMA species codes to NCBI taxon IDs and scientific names.
# Note: OMA species codes are whenever possible identical to UniProt codes.
# Format: OMA code<tab>OMA Taxon ID<tab>NCBI Taxon ID<tab>GTDB genome accession<tab>Scientific name<tab>Genome source<tab>Version/Release
ABSGL	4829	4829	n/a	Absidia glauca	EnsemblGenomes	Ensembl Fungi 48; AG_v1
ACAM1	-644103924	329726	RS_GCF_000018105.1	Acaryochloris marina (strain MBIC 11017)	Genome Reviews	18-MAR-2008 (Rel. 88, Last updated, Version 2)
ACAPL	133434	133434	n/a	Acanthaster planci	Refseq	Refseq; OKI-Apl_1.0; GCF_001949145.1; 08-AUG-2017
ACCPU	-745294054	522306	RS_GCF_000024165.1	Accumulibacter phosphatis (strain UW-1)	EBI	15-MAY-2014 (Rel. 120, Last updated, Version 4)
ACEAZ	-832469519	574087	RS_GCF_000144695.1	Acetohalobium arabaticum (strain ATCC 49924 / DSM 5501 / Z-7288)	EBI	09-JAN-2014 (Rel. 119, Last updated, Version 5)
ACECE	-823481299	720554	RS_GCF_000237085.1	Acetivibrio clariflavus (strain DSM 19732 / NBRC 101661 / EBR45)	EBI	23-JUL-2013 (Rel. 117, Last updated, Version 3)
ACEMN	-709043147	891968	RS_GCF_000266925.1	Acetomicrobium mobile (strain ATCC BAA-54 / DSM 13181 / JCM 12221 / NGA)	EBI	23-JUL-2013 (Rel. 117, Last updated, Version 3)
```
* [`oma-uniprot`](https://omabrowser.org/All/oma-groups.txt.gz) 
```
# Mapping of OMA IDs to Uniprot for OMA release of All.Jun2023
# The mapping is to both TrEML and Swiss-Prot and also to
# any of the protein names and primary and secondary ACs.
# Format: OMA ID<tab>UniProt ID
ARCFU00001	O30234
ARCFU00002	O30233
ARCFU00010	A0A101DYQ4
ARCFU00012	A0A101DFD4
ARCFU00013	A0A101DFE5
ARCFU00014	A0A101DFH3
```