# SynGAP
usage: SynGAP [-h] {dual,triple,custom,genepair,evi,eviplot} ...

| **Module** | **Description** |
| --- | --- |
| _dual_ | Genome annotations polishment for two species |
| _triple_ | Genome annotations polishment for three species |
| _custom_ | Genome annotations polishment for two species within specific synteny blocks |
| _genepair_ | Pairing of high-confidence cross-species homologous genes |
| _evi_ | _EVI_ calculation for each gene pair |
| _eviplot_ | Ranked dotplot for given gene pairs with _EVI_ and highlight specific gene pairs |
| h | Show this help message and exit |

## _dual_
usage: SynGAP dual [-h] --sp1fa SP1FA --sp1gff SP1GFF --sp2fa SP2FA --sp2gff SP2GFF 
           [--sp1 SP1] [--sp2 SP2]
           [--annoType1 ANNOTYPE1] [--annoKey1 ANNOKEY1] [--annoparentKey1 ANNOPARENTKEY1]
           [--annoType2 ANNOTYPE2] [--annoKey2 ANNOKEY2] [--annoparentKey2 ANNOPARENTKEY2]
           [--datatype DATATYPE] [--cscore CSCORE] [--threads THREADS]
           [--process PROCESS] [--evalue EVALUE] [--rank RANK] [--coverage COVERAGE]
           [--kmer1 KMER1] [--kmer2 KMER2] [--outs OUTS] [--intron INTRON]

| **Parameter** | **Description** |
| --- | --- |
| sp1fa | The genome squence file (.fasta format) for species1 [required] |
| sp1gff | The genome annotation file (.gff format) for species1 [required] |
| sp2fa | The genome squence file (.fasta format) for species2 [required] |
| sp2gff | The genome annotation file (.gff format) for species2 [required] |
| sp1 | The short name for species1, e.g. Ath [default: sp1] |
| sp2 | The short name for species2, e.g. Ath [default: sp2] |
| annoType1 | Feature type to extract for species1 [default: mRNA] |
| annoKey1 | Key in the attributes to extract for species1 [default: ID] |
| annoparentKey1 | Parent gene key to group with --primary_only in jcvi [default: Parent] |
| annoType2 | Feature type to extract for species2 [default: mRNA] |
| annoKey2 | Key in the attributes to extract for species2 [default: ID] |
| annoparentKey2 | Parent gene key to group with --primary_only in jcvi [default: Parent] |
| datatype | The type of squences for jcvi, nucl&#124;prot [default: nucl] |
| cscore | C-score cutoff for jcvi [default: 0.7] |
| threads | Number of threads to use [default: 8] |
| process | Process for gapanno, genblastg&#124;miniprot [default: genblastg] |
| evalue | Threshold for evalue in genBlast [default: 1e-5] |
| rank | The number of ranks in genBlast output [default: 5] |
| coverage | Minimum percentage of query gene coverage of the HSP group in the genBlast output [default: 0.5] |
| kmer1 | K-mer size for Indexing in miniprot [default: 5] |
| kmer2 | K-mer size for the second round of chaining in miniprot [default: 4] |
| outs | Threshold of Score for miniprot output [default: 0.95] |
| intron | Max intron size allowed for miniprot output [default: 40k] |
| h | Show this help message and exit |

## _initdb_
usage: SynGAP initdb [-h] --sp SP --file FILE

| **Parameter** | **Description** |
| --- | --- |
| sp | The species type of masterdb to be imported, plant&#124;animal [required] |
| file | The compressed file of masterdb (.tar.gz) to be imported [required] |
| h | Show this help message and exit |

## _master_
usage: SynGAP master [-h] --sp SP
           --sp1fa SP1FA --sp1gff SP1GFF [--sp1 SP1]
           [--annoType1 ANNOTYPE1] [--annoKey1 ANNOKEY1] [--annoparentKey1 ANNOPARENTKEY1]
           [--datatype DATATYPE] [--cscore CSCORE] [--threads THREADS]
           [--process PROCESS] [--evalue EVALUE] [--rank RANK] [--coverage COVERAGE]
           [--kmer1 KMER1] [--kmer2 KMER2] [--outs OUTS] [--intron INTRON]

| **Parameter** | **Description** |
| --- | --- |
| sp | The species type of the polished object, plant&#124;animal [required] |
| sp1fa | The genome squence file (.fasta format) for species1 [required] |
| sp1gff | The genome annotation file (.gff format) for species1 [required] |
| sp1 | The short name for species1, e.g. Ath [default: sp1] |
| annoType1 | Feature type to extract for species1 [default: mRNA] |
| annoKey1 | Key in the attributes to extract for species1 [default: ID] |
| annoparentKey1 | Parent gene key to group with --primary_only in jcvi [default: Parent] |
| datatype | The type of squences for jcvi, nucl&#124;prot [default: nucl] |
| cscore | C-score cutoff for jcvi [default: 0.7] |
| threads | Number of threads to use [default: 8] |
| process | Process for gapanno, genblastg&#124;miniprot [default: genblastg] |
| evalue | Threshold for evalue in genBlast [default: 1e-5] |
| rank | The number of ranks in genBlast output [default: 5] |
| coverage | Minimum percentage of query gene coverage of the HSP group in the genBlast output [default: 0.5] |
| kmer1 | K-mer size for Indexing in miniprot [default: 5] |
| kmer2 | K-mer size for the second round of chaining in miniprot [default: 4] |
| outs | Threshold of Score for miniprot output [default: 0.95] |
| intron | Max intron size allowed for miniprot output [default: 40k] |
| h | Show this help message and exit |

## _triple_
usage: SynGAP triple [-h] --sp1fa SP1FA --sp1gff SP1GFF --sp2fa SP2FA --sp2gff SP2GFF 
           --sp3fa SP3FA --sp3gff SP3GFF [--sp1 SP1] [--sp2 SP2] [--sp3 SP3]
           [--annoType1 ANNOTYPE1] [--annoKey1 ANNOKEY1] [--annoparentKey1 ANNOPARENTKEY1]
           [--annoType2 ANNOTYPE2] [--annoKey2 ANNOKEY2] [--annoparentKey2 ANNOPARENTKEY2]
           [--annoType3 ANNOTYPE3] [--annoKey3 ANNOKEY3] [--annoparentKey3 ANNOPARENTKEY3]
           [--datatype DATATYPE] [--cscore CSCORE] [--threads THREADS]
           [--process PROCESS] [--evalue EVALUE] [--rank RANK] [--coverage COVERAGE]
           [--kmer1 KMER1] [--kmer2 KMER2] [--outs OUTS] [--intron INTRON]

| **Parameter** | **Description** |
| --- | --- |
| sp1fa | The genome squence file (.fasta format) for species1 [required] |
| sp1gff | The genome annotation file (.gff format) for species1 [required] |
| sp2fa | The genome squence file (.fasta format) for species2 [required] |
| sp2gff | The genome annotation file (.gff format) for species2 [required] |
| sp3fa | The genome squence file (.fasta format) for species3 [required] |
| sp3gff | The genome annotation file (.gff format) for species3 [required] |
| sp1 | The short name for species1, e.g. Ath [default: sp1] |
| sp2 | The short name for species2, e.g. Ath [default: sp2] |
| sp3 | The short name for species3, e.g. Ath [default: sp3] |
| annoType1 | Feature type to extract for species1 [default: mRNA] |
| annoKey1 | Key in the attributes to extract for species1 [default: ID] |
| annoparentKey1 | Parent gene key to group with --primary_only in jcvi [default: Parent] |
| annoType2 | Feature type to extract for species2 [default: mRNA] |
| annoKey2 | Key in the attributes to extract for species2 [default: ID] |
| annoparentKey2 | Parent gene key to group with --primary_only in jcvi [default: Parent] |
| annoType3 | Feature type to extract for species3 [default: mRNA] |
| annoKey3 | Key in the attributes to extract for species3 [default: ID] |
| annoparentKey3 | Parent gene key to group with --primary_only in jcvi [default: Parent] |
| datatype | The type of squences for jcvi, nucl&#124;prot [default: nucl] |
| cscore | C-score cutoff for jcvi [default: 0.7] |
| threads | Number of threads to use [default: 8] |
| process | Process for gapanno, genblastg&#124;miniprot [default: genblastg] |
| evalue | Threshold for evalue in genBlast [default: 1e-5] |
| rank | The number of ranks in genBlast output [default: 5] |
| coverage | Minimum percentage of query gene coverage of the HSP group in the genBlast output [default: 0.5] |
| kmer1 | K-mer size for Indexing in miniprot [default: 5] |
| kmer2 | K-mer size for the second round of chaining in miniprot [default: 4] |
| outs | Threshold of Score for miniprot output [default: 0.95] |
| intron | Max intron size allowed for miniprot output [default: 40k] |
| h | Show this help message and exit |

## _custom_
usage: SynGAP custom [-h] --sp1fa SP1FA --sp1gff SP1GFF --sp2fa SP2FA --sp2gff SP2GFF
           --custom_anchors CUSTOM_ANCHORS
           [--sp1 SP1] [--sp2 SP2]
           [--annoType1 ANNOTYPE1] [--annoKey1 ANNOKEY1] [--annoparentKey1 ANNOPARENTKEY1]
           [--annoType2 ANNOTYPE2] [--annoKey2 ANNOKEY2] [--annoparentKey2 ANNOPARENTKEY2]
           [--datatype DATATYPE] [--cscore CSCORE] [--threads THREADS]
           [--process PROCESS] [--evalue EVALUE] [--rank RANK] [--coverage COVERAGE]
           [--kmer1 KMER1] [--kmer2 KMER2] [--outs OUTS] [--intron INTRON]

| **Parameter** | **Description** |
| --- | --- |
| sp1fa | The genome squence file (.fasta format) for species1 [required] |
| sp1gff | The genome annotation file (.gff format) for species1 [required] |
| sp2fa | The genome squence file (.fasta format) for species2 [required] |
| sp2gff | The genome annotation file (.gff format) for species2 [required] |
| sp1 | The short name for species1, e.g. Ath [default: sp1] |
| sp2 | The short name for species2, e.g. Ath [default: sp2] |
| custom_anchors | Choose the self-defined Syntenic anchors file [required] |
| annoType1 | Feature type to extract for species1 [default: mRNA] |
| annoKey1 | Key in the attributes to extract for species1 [default: ID] |
| annoparentKey1 | Parent gene key to group with --primary_only in jcvi [default: Parent] |
| annoType2 | Feature type to extract for species2 [default: mRNA] |
| annoKey2 | Key in the attributes to extract for species2 [default: ID] |
| annoparentKey2 | Parent gene key to group with --primary_only in jcvi [default: Parent] |
| datatype | The type of squences for jcvi, nucl&#124;prot [default: nucl] |
| cscore | C-score cutoff for jcvi [default: 0.7] |
| threads | Number of threads to use [default: 8] |
| process | Process for gapanno, genblastg&#124;miniprot [default: genblastg] |
| evalue | Threshold for evalue in genBlast [default: 1e-5] |
| rank | The number of ranks in genBlast output [default: 5] |
| coverage | Minimum percentage of query gene coverage of the HSP group in the genBlast output [default: 0.5] |
| kmer1 | K-mer size for Indexing in miniprot [default: 5] |
| kmer2 | K-mer size for the second round of chaining in miniprot [default: 4] |
| outs | Threshold of Score for miniprot output [default: 0.95] |
| intron | Max intron size allowed for miniprot output [default: 40k] |
| h | Show this help message and exit |

## _genepair_
usage: SynGAP genepair [-h] --sp1fa SP1FA --sp1gff SP1GFF --sp2fa SP2FA --sp2gff SP2GFF 
           [--sp1 SP1] [--sp2 SP2]
           [--annoType1 ANNOTYPE1] [--annoKey1 ANNOKEY1] [--annoparentKey1 ANNOPARENTKEY1]
           [--annoType2 ANNOTYPE2] [--annoKey2 ANNOKEY2] [--annoparentKey2 ANNOPARENTKEY2]
           [--datatype DATATYPE] [--cscore CSCORE] [--evalue EVALUE] [--iTAK ITAK] [--threads THREADS]

| **Parameter** | **Description** |
| --- | --- |
| sp1fa | The genome squence file (.fasta format) for species1 [required] |
| sp1gff | The genome annotation file (.gff format) for species1 [required] |
| sp2fa | The genome squence file (.fasta format) for species2 [required] |
| sp2gff | The genome annotation file (.gff format) for species2 [required] |
| sp1 | The short name for species1, e.g. Ath [default: sp1] |
| sp2 | The short name for species2, e.g. Ath [default: sp2] |
| annoType1 | Feature type to extract for species1 [default: mRNA] |
| annoKey1 | Key in the attributes to extract for species1 [default: ID] |
| annoparentKey1 | Parent gene key to group with --primary_only in jcvi [default: Parent] |
| annoType2 | Feature type to extract for species2 [default: mRNA] |
| annoKey2 | Key in the attributes to extract for species2 [default: ID] |
| annoparentKey2 | Parent gene key to group with --primary_only in jcvi [default: Parent] |
| datatype | The type of squences for jcvi, nucl&#124;prot [default: nucl] |
| cscore | C-score cutoff for jcvi [default: 0.7] |
| evalue | Threshold for evalue in two-way blast [default: 1e-2] |
| iTAK | Perform iTAK to identify TFs and kinases (only for plants), yse&#124;no [default: no] |
| threads | Number of threads to use [default: 8] |
| h | Show this help message and exit |

## _evi_
usage: SynGAP evi [-h] --genepair GENEPAIR --sp1exp SP1EXP --sp2exp SP2EXP
           [--weight WEIGHT] [--format FORMAT]

| **Parameter** | **Description** |
| --- | --- |
| genepair | The genepair file (tab-divided) for _EVI_ counting [required] |
| sp1exp | The expression file (tab-divided) for species1 [required] |
| sp2exp | The expression file (tab-divided) for species2 [required] |
| weight | The weight of three indexes in _EVI_ calulation (ML:FC:PCC) [default=1:1:4] |
| format | The format of output figure |
| h | Show this help message and exit |

## _eviplot_
usage: SynGAP eviplot [-h] --EVI EVI [--highlightid HIGHLIGHTID] [--highlightcolor HIGHLIGHTCOLOR]
           --outgraph OUTGRAPH [--figsize FIGSIZE] [--fontsize FONTSIZE]

| **Parameter** | **Description** |
| --- | --- |
| EVI | The _EVI_ file (tab-divided) for ploting [required] |
| highlightid | The id list file (tab-divided) for highlightid |
| highlightcolor | The color for highlight label [default=red] |
| outgraph | The output graph file (output format determined by the file suffix) [required] |
| figsize | The size of output graph (LengthxWidth) [default=10x5] |
| fontsize | The font size of output graph [default=10] |
| h | Show this help message and exit |


