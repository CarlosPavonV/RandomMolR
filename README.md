# RandomMolR
Random scripts for molecular data edition in R:<br/><br/>
BPPformat creates a sequence file in BPP format based on a concatenated FASTA alignment and a partitions file in raxml format.<br/><br/>
collapse_conStruct is used to edit the input files used by conStruct (https://github.com/gbradburd/conStruct) by collapsing localities that are geographically close to each other based on a distance threshold and updating the allele frequency, coordinate, and geographic distance matrices accordingly.<br/><br/>
CreateIndLoc splits a concatenated FASTA alignment into individual locus alignments based on a partitions file in raxml format.<br/><br/>
DartPart takes a genlight object (such as that used in dartR: https://github.com/green-striped-gecko/dartR) and creates a partitions file in raxml or nexus format. Optionally, it writes files with the individual alignments of each locus.<br/><br/>
GBclean shortens the names of sequences downloaded from GenBank in a somewhat flexible way.<br/><br/>
seattleR is an R version of SeATTLE (Sequence Alignment Transformation into a Table for Later Edition), which takes a concatenated FASTA alignment and a partitions file to give you a data frame in which columns correspond to loci. Among other things, this table can be used as input for phrynomics (https://github.com/bbanbury/phrynomics/) to extract and manipulate SNPs.<br/><br/>
TerGaps2N replaces terminal gaps with a symbol denoting missing data (e.g. N or ?) in a FASTA alignemnt. This is useful when downstream analysis involves software that distinguishes between gaps and missing data. It also makes it easier to look for undesired internal gaps.<br/><br/>
