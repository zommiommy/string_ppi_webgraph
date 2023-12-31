# string_ppi_webgraph
Experiments on merging and compressing the whole inter-species string ppi graph using webgraph-rs


#### MANDATORY
This will create many big files, please run this to increase the limit:
```bash
sudo sysctl -w fs.file-max=100000000
sudo ulimit -sn 100000
```
and set a tmp directory with >100GB of free space like:
```bash
export TMPDIR="/dfd/tmp"
```
moreover, all these operaion might require time, so you can enable cpu specific
optimizations like:
```bash
export RUSTFLAGS="-C target-cpu=native"
```

### How to
**To build the vocabulary run:**
```bash
cargo run --release --bin vocab_builder
```
this will create a `vocab.tsv` file, where the node names of each node id can 
be found at the line with the corresponding index, and `vocab.sorted.tsv` which
is a list of `node_name,node_id`.

**Then, to build the graph run:**
```bash
cargo run --release --bin graph_builder
```
This will create `res.graph` which is the compressed graph, and `res.properties`
which contains metadata needed to read the graph and some other stats.

**To convert a webgraph to tsv**, inside the [webgraph-rs](https://github.com/vigna/webgraph-rs) repository, use:
```bash
cargo run --release --bin to_csv /path/to/webgraph/
```
be careful that you must not include the `.graph` suffix.

*Currently there is a bug in the paralle compression so we will add `-j1` in the
following commands, once the bug is fixed it can be removed.*

**To run layered label propagation**, inside the [webgraph-rs](https://github.com/vigna/webgraph-rs) repository, use:
```bash
# make the graph undirected and without selfoops, this will create res.simple.graph and res.simple.properties.
cargo run --release --bin simplify res -j1 --batch-size=100000000
# build elias-fano for random access, this will create res.simple.ef
cargo run --release --bin build_eliasfano res.simple
# run the llp this will create res.simple.llp, and 12 label files which 
# can be used for comunity detection
cargo run --release --bin llp res.simple
# apply the llp permutation to the original (possibly directed) graph.
# this will create the, hopefully, better compressed graph `res.comp.graph` and
# `res.comp.properties`
cargo run --release --bin perm res res.comp res.simple.llp -j1 --batch-size=100000000
```
`batch-size` is how many arcs (16 bytes) we will keep in memory to sort and then
dump on a file for the external sort. The bigger the value, the faster the 
compression will be, but it will will use more RAM.

# Sources

```python
urls = [
    "https://omabrowser.org/All/oma-uniprot.txt.gz",
    "https://omabrowser.org/All/oma-species.txt",
    "https://omabrowser.org/All/oma-groups.txt.gz",

    "http://eggnog6.embl.de/download/eggnog_6.0/e6.og2seqs_and_species.tsv",

    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz",

    "https://kg-hub.berkeleybop.io/kg-obo/ncbitaxon/2023-09-14/ncbitaxon_kgx_tsv.tar.gz",
    "https://kg-hub.berkeleybop.io/kg-obo/go/2023-03-06/go_kgx_tsv.tar.gz",
    "https://kg-hub.berkeleybop.io/kg-obo/bto/2021-10-26/bto_kgx_tsv.tar.gz",
    "https://kg-hub.berkeleybop.io/kg-obo/cl/2023-10-19/cl_kgx_tsv.tar.gz",
    "https://kg-hub.berkeleybop.io/kg-obo/doid/2023-10-21/doid_kgx_tsv.tar.gz",
    "https://kg-hub.berkeleybop.io/kg-obo/fbcv/2023-10-19/fbcv_kgx_tsv.tar.gz",
    "https://kg-hub.berkeleybop.io/kg-obo/hp/2023-10-09/hp_kgx_tsv.tar.gz",
    "https://kg-hub.berkeleybop.io/kg-obo/mp/2021-11-04/mp_kgx_tsv.tar.gz",
    "https://kg-hub.berkeleybop.io/kg-obo/zp/2022-12-17/zp_kgx_tsv.tar.gz",

    "https://stringdb-downloads.org/download/protein.info.v12.0.txt.gz",
    "https://stringdb-downloads.org/download/protein.aliases.v12.0.txt.gz",
    "https://stringdb-downloads.org/download/protein.enrichment.terms.v12.0.txt.gz",
    "https://stringdb-downloads.org/download/protein.links.full.v12.0.txt.gz",
]

from downloaders import BaseDownloader
downloader = BaseDownloader()
for url in urls:
    _ = downloader.download(url)
```

Needed files from [string ppi](https://string-db.org/cgi/download):
* [`protein.info.v12.0.txt.gz`](https://stringdb-downloads.org/download/protein.info.v12.0.txt.gz)
```
#string_protein_id      preferred_name  protein_size    annotation                                                                   
23.BEL05_00025  BEL05_00025     429     Hypothetical protein; Incomplete; too short partial abutting assembly gap; missing start and stop; Derived by automated computational analysis using gene prediction method: GeneMarkS+.                                          
23.BEL05_00030  OEG73659.1      203     2-keto-4-pentenoate hydratase; Derived by automated computational analysis using gene prediction method: Protein Homology.                                                                                                        
23.BEL05_00035  OEG73660.1      1140    Hybrid sensor histidine kinase/response regulator; Derived by automated computational analysis using gene prediction method: Protein Homology.                                                                                    
23.BEL05_00040  acsA    650     acetate--CoA ligase; Catalyzes the conversion of acetate into acetyl-CoA (AcCoA), an essential intermediate at the junction of anabolic and catabolic pathways. AcsA undergoes a two-step reaction. In the first half reaction, AcsA combines acetate with ATP to form acetyl-adenylate (AcAMP) intermediate. In the second half reaction, it can then transfer the acetyl group from AcAMP to the sulfhydryl group of CoA, forming the product AcCoA; Belongs to the ATP-dependent AMP-binding enzyme family.      
23.BEL05_00045  OEG73662.1      579     ATP-dependent helicase; Derived by automated computational analysis using gene prediction method: Protein Homology.                                                                                                               
23.BEL05_00050  OEG73663.1      85      Glutaredoxin; Functions as an electron carrier in the glutathione-dependent synthesis of deoxyribonucleotides by the enzyme ribonucleotide reductase; also involved in reducing some disulfides in a coupled system with glutathione reductase; Derived by automated computational analysis using gene prediction method: Protein Homology.                            
23.BEL05_00055  OEG73664.1      380     Porin; Derived by automated computational analysis using gene prediction method: Protein Homology.                                                                                                                                
23.BEL05_00070  cpoB    241     Tol-pal system protein YbgF; Mediates coordination of peptidoglycan synthesis and outer membrane constriction during cell division; Belongs to the CpoB family.                                                                           
23.BEL05_00075  pal     180     Peptidoglycan-associated lipoprotein; Part of the Tol-Pal system, which plays a role in outer membrane invagination during cell division and is important for maintaining outer membrane integrity.                                       
```
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
* [`protein.aliases.v12.0.txt.gz`](https://stringdb-downloads.org/download/protein.aliases.v12.0.txt.gz)
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
* [`e6.og2seqs_and_species.tsv`](http://eggnog6.embl.de/download/eggnog_6.0/e6.og2seqs_and_species.tsv) The first column is the taxonomy ID which is the last common ancestor, which is the NCBI taxon. The second column is the name of the ortholog group. third and fourth column is the number of species and proteins respectively in the group. fifth column is the ids of the species involved in the group. and the sixth column is the STRING ids of the proteins (these are more orthologs).
```
314146  4R1PH   3       3       591936,61622,336983     336983.ENSCANP00000003193,591936.ENSPTEP00000002293,61622.ENSRROP00000006056
314146  4R1PI   2       2       61622,336983    336983.ENSCANP00000035642,61622.ENSRROP00000039496
314146  4R1PJ   2       4       61622,336983    336983.ENSCANP00000009298,336983.ENSCANP00000014889,61622.ENSRROP00000001152,61622.ENSRROP00000003608
314146  4R1PK   3       3       591936,61622,61621      591936.ENSPTEP00000029268,61621.ENSRBIP00000000020,61622.ENSRROP00000043585
314146  4R1PM   2       2       591936,61621    591936.ENSPTEP00000024160,61621.ENSRBIP00000003159
314146  4R1PN   3       4       61621,61622,336983      336983.ENSCANP00000007271,336983.ENSCANP00000011724,61621.ENSRBIP00000019800,61622.ENSRROP00000015263
314146  4R1PP   2       2       61622,336983    336983.ENSCANP00000024027,61622.ENSRROP00000027054
314146  4R1PQ   2       2       61622,336983    336983.ENSCANP00000020242,61622.ENSRROP00000005057
314146  4R1PR   2       2       591936,336983   336983.ENSCANP00000031012,591936.ENSPTEP00000036989
314146  4R1PS   2       2       61621,336983    336983.ENSCANP00000038913,61621.ENSRBIP00000035373
```
* [`uniprot_sprot.dat.gz`](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz)
```
ID   001R_FRG3G              Reviewed;         256 AA.
AC   Q6GZX4;
DT   28-JUN-2011, integrated into UniProtKB/Swiss-Prot.
DT   19-JUL-2004, sequence version 1.
DT   08-NOV-2023, entry version 44.
DE   RecName: Full=Putative transcription factor 001R;
GN   ORFNames=FV3-001R;
OS   Frog virus 3 (isolate Goorha) (FV-3).
OC   Viruses; Varidnaviria; Bamfordvirae; Nucleocytoviricota; Megaviricetes;
OC   Pimascovirales; Iridoviridae; Alphairidovirinae; Ranavirus; Frog virus 3.
```
* [`go_kgx_tsv.tar.gz`](https://kg-hub.berkeleybop.io/kg-obo/go/2023-03-06/go_kgx_tsv.tar.gz)
`go_kgx_tsv_edges.tsv`
`go_kgx_tsv_nodes.tsv`

* [`ncbitaxon_kgx_tsv.tar.gz`](https://kg-hub.berkeleybop.io/kg-obo/ncbitaxon/2023-09-14/ncbitaxon_kgx_tsv.tar.gz)
`ncbitaxon_kgx_tsv_edges.tsv`
```
id      subject predicate       object  relation        knowledge_source                                                                  
urn:uuid:8ba9a91f-c534-4371-a789-8a3e676ca1b7   NCBITaxon:1252590       biolink:subclass_of     NCBITaxon:173037        rdfs:subClassOf NCBITAXON 2023-09-14                                                                                                                        
urn:uuid:4e7b5a22-12ba-4d0d-bca1-8d60b583b02e   NCBITaxon:1252591       biolink:subclass_of     NCBITaxon:173037        rdfs:subClassOf NCBITAXON 2023-09-14                                                                                                                        
urn:uuid:62b750fc-4633-4860-8e5e-cb0456f7a71c   NCBITaxon:1252592       biolink:subclass_of     NCBITaxon:173037        rdfs:subClassOf NCBITAXON 2023-09-14                                                                                                                        
urn:uuid:7859d3ce-b49b-44c6-9160-7b5af01f5dc8   NCBITaxon:1252593       biolink:subclass_of     NCBITaxon:2633401       rdfs:subClassOf NCBITAXON 2023-09-14                                                                                                                        
urn:uuid:0f86d260-c6e3-41d8-9b38-98733e23cedc   NCBITaxon:1252594       biolink:subclass_of     NCBITaxon:173037        rdfs:subClassOf NCBITAXON 2023-09-14                                                                                                                        
urn:uuid:abeae0d5-72a5-43de-950d-2f11e3e15c93   NCBITaxon:1252595       biolink:subclass_of     NCBITaxon:632588        rdfs:subClassOf NCBITAXON 2023-09-14                                                                                                                        
urn:uuid:6f751bf1-3de7-4df3-a3f2-0a4576e821c1   NCBITaxon:1252596       biolink:subclass_of     NCBITaxon:632221        rdfs:subClassOf NCBITAXON 2023-09-14                                                                                                                        
urn:uuid:27b0201b-9d86-4d55-908e-7defcbb34da5   NCBITaxon:1252597       biolink:subclass_of     NCBITaxon:632256        rdfs:subClassOf NCBITAXON 2023-09-14                                                                                                                        
urn:uuid:270736fb-298e-4601-9513-2615b6bbc368   NCBITaxon:1252598       biolink:subclass_of     NCBITaxon:274104        rdfs:subClassOf NCBITAXON 2023-09-14         
```
`ncbitaxon_kgx_tsv_nodes.tsv`
```
id      category        name    description     xref    provided_by     synonym iri     same_as                                           
NCBITaxon:1252590       biolink:OrganismalEntity        Ichneumonidae sp. BOLD:AAG9159          GC_ID:1 ncbitaxon.json          http://purl.obolibrary.org/obo/NCBITaxon_1252590                                                                                                    
NCBITaxon:1252591       biolink:OrganismalEntity        Ichneumonidae sp. ASWA193-08            GC_ID:1 ncbitaxon.json          http://purl.obolibrary.org/obo/NCBITaxon_1252591                                                                                                    
NCBITaxon:1252592       biolink:OrganismalEntity        Ichneumonidae sp. BOLD:AAH1531          GC_ID:1 ncbitaxon.json          http://purl.obolibrary.org/obo/NCBITaxon_1252592                                                                                                    
NCBITaxon:1252593       biolink:OrganismalEntity        Stibeutes sp. BOLD:ABZ5036              GC_ID:1 ncbitaxon.json          http://purl.obolibrary.org/obo/NCBITaxon_1252593                                                                                                    
NCBITaxon:1252594       biolink:OrganismalEntity        Ichneumonidae sp. BOLD:AAH2156          GC_ID:1 ncbitaxon.json          http://purl.obolibrary.org/obo/NCBITaxon_1252594                                                                                                    
NCBITaxon:1252595       biolink:OrganismalEntity        Mesochorinae sp. BOLD:AAH1799           GC_ID:1 ncbitaxon.json          http://purl.obolibrary.org/obo/NCBITaxon_1252595                                                                                                    
NCBITaxon:1252596       biolink:OrganismalEntity        Alysiinae sp. ASWAT849-08               GC_ID:1 ncbitaxon.json          http://purl.obolibrary.org/obo/NCBITaxon_1252596                                                                                                    
NCBITaxon:1252597       biolink:OrganismalEntity        Campopleginae sp. BOLD:AAF8320          GC_ID:1 ncbitaxon.json          http://purl.obolibrary.org/obo/NCBITaxon_1252597                                                                                                    
NCBITaxon:1252598       biolink:OrganismalEntity        Cryptinae sp. BOLD:ABZ0478              GC_ID:1 ncbitaxon.json          http://purl.obolibrary.org/obo/NCBITaxon_1252598```     