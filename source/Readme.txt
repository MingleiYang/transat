* Making Transat from source:
	Transat requires the alignment shuffling perl script from RNAz. RNAz can be found at http://www.tbi.univie.ac.at/~wash/RNAz/
	
	The makefile is located in the Transat directory. However, you must first edit the makefile to specify where the rnazRandomizeAln.pl script can be found on your machine. Set value the RNAZ_SHUFFLER_LOC variable in the makefile to the appropriate location in your filesystem.
	
	Type 'make all' in the Transat directory to compile Transat
		
* Transat generates a data table with an entry for every helix in the alignment
	* notes:
	All helices in this table are unique (unlike previous tables)
	P-values are calculated by shuffling entire alignment and finding all helices.
	P-values are raw (not B-H corrected)
	Order of the columns is very different from previous tables. Hopefully the column descriptions are relatively self-explanatory.
	
Transat/Transat -fasta <fasta alignment file> -ct <consensus structure .ct file> -tree <tree in newick format>
or:
Transat/Transat -fasta <fasta alignment file, with last entry as the structure in dot-bracket notation> -tree <tree in newick format>
options:

-fasta [filename] : file containing alignment, in aligned fasta format
-ct [filename] : file containing known structure. If a ct file not specified, then Transat will assume that the last line of the alignment file contains the structure in dot-bracket notation.
-tree [filename] : file containing phylogenetic tree in newick format.
-minSL [int] : sets the minimum helix length. Only competing helices with a length greater than this value will be stored (default: 8, but for all my rfam analyses, i've been using 3)
-randomize [int] : sets the number of shuffled alignments are used in the null distributions (default = 1000)
-realign : if this flag is present, Transat does a realignment step before producing a shuffled alignments. Uses t_coffee for realignment step (t_coffee must be on PATH!)
-noDB : if this flag is present, Transat will not print the dot-bracket structure of each helix to the output file.
-noStruct : use this flag if you would like to run Transat without specifying a known structure for the alignment. -ct overrides this flag.

see example directory for more information...
