/*
 * Copyright: Nicholas P Wiebe (2008-2009) and Irmtraud M Meyer (2008-2009)
 * License: licensed under the GNU General Public License version 3 (GPLv3)
 * Alignment.h
 *
 * Class for storing an alignment, and associated information
 * Stores sequences in as aligned and unaligned (gapless) strings
 * Also stores the consensus structure for the alignment, and labels
 * the consensus helices of the structure.
 *
 * The method alignmentGStats finds and stores all competing helices in the sequences
 * of the alignment. Many of the other methods in this class require it to be run first.
 * TODO: offload competing helix finding onto HelixFinder Class
 */

#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_

#include <vector>
#include <string>
#include <set>
#include "StatsWrapper.h"
#include "CompetingHelix.h"
#include "InterestingRegion.h"
#include "HelixGroup.h"
#include "HelixCore.h"
#include "Helix.h"
#include "UTMatrix.h"
#include "HelixFinder.h"

//Some of these are unused, I think...
#define MAX_SEQ_L  4500    // maximal sequence length
#define MAX_WORD_L  100    // maximal length of a char* word
#define MAX_LINE_L 1000    // maximal length of a line
#define MIN_DIST      3    // minimal sequence distance between two base-pairing letters

using namespace std;

//forward declaration... Maybe this is avoidable?
class Tree;

class Alignment {
public:
	/*
	 * For this constructor, the structure must be present in dot bracket notation as the last entry
	 * of the fasta alignment file
	 *
	 */
	Alignment(string& filename);

	//alignment file is an aligned fasta file
	//struct file is a sequence file
	Alignment(string& filename, string& structFilename);

	//create an alignment from vectors specifying sequence names, aligned sequences, and consensus structure
	Alignment(vector<string> & names, vector<string> & alignment, vector<int> & structure);

	//blank constructor - do nothing
	Alignment();

	virtual ~Alignment();

	/**
	 * Finds all competing helices from each sequence in the alignment and stores them
	 * in competingHelices.
	 * Also stores g-stats for each competing-helix/true helix pair in StatsMatrix
	 */
	void alignmentGStats(const double eStack[4][4][4][4]);
	/**
	 * Finds all competing helices from a specified sequence in the alignment and stores them
	 * in competingHelices.
	 * Also stores g-stats for each competing-helix/true helix pair in StatsMatrix
	 */
	void sequenceGStats(const double eStack[4][4][4][4], const int seqIndex);

	/**
	 * various methods for printing information about the alignment, mostly for debugging
	 * purposes.
	 */
	string alignmentString();
	string seqString();
	string alignedStructString();
	string structString();
	string alignedlabelString();
	string labelString();
	string competingHelicesString();

	/**
	 * prints the structure of a single consensus helix in dot-bracket notation
	 * deprecated in favor of TrueHelixString (though i suspect they both work
	 * exactly the same).
	 */
	string printTrueHelix(int i);



	/* added by rgoya */
	/** deprecated */
	string gappedStructString(); // similar to structString but shifts according to gaps,
								// makes visualizing the helices much easier on the eyes
	/** deprecated */
	string competingHelicesHistogram(bool naive); // calculates values for the histogram of competing helices hits
												  // if(naive) -> mark both true and competing helix parts
												  // else      -> only mark competing helix part
	/** deprecated */
	string competingHelicesHistogramTESTING(bool naive);
	/** deprecated */
	string competingHelicesHistogramPerHelixTESTING(int trueHelix);
	/** deprecated */
	vector<vector<vector<double*>* >* > competingHelicesHistogramPerHelixVector(int trueHelix);
	/*
	 * Returns the number of sequences we have to work with
	 * (mostly I just use alignedSeqs.size() instead... -nick)
	 */
	unsigned int numSeqs();

	/*
	 * Returns the number of true helices in the structure
	 */
	unsigned int numTrueHelices();

	/*
	 * Returns the number of competing helices affecting
	 * a specific true helix given as a parameter
	 */
	unsigned int numCompetingOnTrueHelix(int trueHelix);

	/*
	 * Returns the number of competing helices affecting
	 * a specific true helix in a specific sequence given
	 * as parameters
	 */
	unsigned int numCompetingOnSeqTrueHelix(int seqIndex, int trueHelix);

	string printCompetingHelices(int trueHelix);
	/* end added by rgoya */

	/*
	 * computes and prints to standard out the information
	 * content of each column in the alignment
	 * (a measure of primary sequence conservation)
	 * high = highly conserved
	 *
	 * deprecated
	 */
	void alignmentInformationContent();

	/*
	 * A more straightforward measure of sequence conservation
	 *
	 * deprecated
	 */
	void sequenceConservation();

	/*
	 * sequence conservation for alignment position i
	 *
	 * deprecated
	 */
	double seqConservationAtPosition(int i);

	/*
	 * Prints to std out the number of competing helices for each sequence
	 * (for making boxplot)
	 *
	 * deprecated
	 */
	void printCompetingHelixCount();

	/*
	 * Prints Cis values for each true helix, for every sequence.
	 *
	 * deprecated
	 */
	void printCisTransGValuesPerHelix();


	// Probable repeated functions?
	/*
	 * Prints Cis and Trans Values for each sequence
	 *
	 * deprecated
	 */
	void printCisTransGValuesPerSequence();

	/*
	 * Print Cis and Trans values (over the whole seqeunce) for each sequence in the alignment
	 *
	 * deprecated
	 */
	void printCisTransValues();


	/*
	 * deprecated
	 */
	void findInterestingRegions(int trueHelix);
	/*
	 * deprecated
	 */
	void rankInterestingRegions();

	/*
	 * deprecated
	 */
	void printCompetingHelixMidpoints();

	/**
	 * Constructs 'Cores' of competing helices, i.e. regions in the alignment
	 * where several sequences have competing helices.
	 *
	 * deprecated
	 */
	void makeCompHelixCores();

	/**
	 * deprecated
	 */
	void printCoreLogLikelihoods(Tree & treeRoot);

	/**
	 * deprecated
	 */
	void printTrueHelixLogLikelihoods(Tree & treeRoot);

	/**
	 * deprecated
	 */
	void printCompetingHelixLogLikelihoods(Tree & treeRoot);

	/**
	 * prints the midpoints for each competing helix and true helix for each competing helix/true helix
	 * interaction pair, and the competing helix log-ratio
	 *
	 * deprecated
	 */
	void printCompetingHelixTrueHelixMidpoints(Tree & treeRoot, const double eStack[4][4][4][4]);

	/**
	 * Shuffles columns of the alignment.
	 * Only shuffles columns that are not part of true helices.
	 *
	 * Standard non-fancy shuffling procedure. This may not be optimal...
	 */
	void shuffleAlignment();

	/**
	 * Calculates p values for the each competing helix.
	 * Simulates the null distribution of competing helices on individual
	 * true helices by shuffling the alignment.
	 *
	 * randomizedTrials := number of randomized alignments to generate for each
	 * true helix
	 *
	 * deprecated in favor of calculatePValues(Tree, string, unsigned int) below
	 */
	void calculatePValues(Tree & tree, unsigned int randomizedTrials);

	/**
	 * Calculates p values for the each competing helix.
	 * Simulates the null distribution of competing helices on individual
	 * true helices by realigning and shuffling the alignment.
	 *
	 * randomizedTrials := number of randomized alignments to generate for each
	 * true helix
	 */
	void calculatePValues(Tree & tree, string & treeFile, unsigned int randomizedTrials);

	/**
	 * deprecated... used for debugging
	 */
	void nickTemp(Tree & tree, unsigned int randomizedTrials);


	/**
	 * Calculates the percentage of competing helices from that participate in core regions
	 *
	 * deprecated
	 */
	double calculatePercentOfCompetingHelicesInCores();

	/**
	 * deprecated
	 */
	double meanValidBPs(CompetingHelix & h, int seqIndex);

	/**
	 * Prints the alignment region of a given competing helix
	 */
	string printAlignmentAtCompetingHelix(CompetingHelix & h, int seqIndex);

	/**
	 * Print data table on consensus helices to stdout
	 */
	void printTrueHelixStats(Tree & t);

	/**
	 *
	 */
	string TrueHelixString(int th_index);

	/**
	 * gets the alignment in Fasta format
	 * @param aligned if true, return the gapped sequences (FastaAln format)
	 *
	 */
	string FastaFormat(bool aligned = true) const;

	/**
	 * gets the alignment in ClustalW format
	 */
	string ClustalWFormat();

	/**
	 * gets the specified columns in ClustalW format
	 */
	string ClustalWFormat(const vector<int> & columns);

	/**
	 * returns fraction of canonical bp for basepair at pos5-pos3
	 * TODO: rename to canonicalBPPercent
	 */
	double consensusBPPercent(int pos5, int pos3);

	/**
	 * returns covariance statistic for basepair at pos5-pos3
	 */
	double getCovariance(int pos5, int pos3);

	/**
	 * returns the sequence conservation at specified alignment position
	 * as measured by mean percent identity
	 */
	double getSeqCons(unsigned int pos);


	/**
	 * Calculates the logLikelihood for a single pair of (alignment) positions
	 * If this has been calculated before, then this is a simple table lookup
	 * If not, then it calculates likelihoods under the paired and singleton models
	 * of evolution using the felsenstein algorithm (O(N), where N is the number of sequences
	 * in the alignment).
	 */
	double logLikelihood(int pos5, int pos3, Tree & tree);

	double logPairedLikelihood(int pos5, int pos3, Tree & tree);
	double logUnpairedLikelihood(int pos, Tree & tree);


	/**
	 * re-aligns the sequences in the alignment within the interval [begin, end)
	 * (Present in Alignment class instead of ShuffledAlignment to help with testing)
	 */
	vector<string> realignInterval(int begin, int end, string & treeFilename) const;

	vector<string> realignTcoffee(Tree & tree) const;
	/**
	 * prints a matrix containing log likelihood ratios for every pair of positions
	 * in the alignment, suitable for importing into R to make an image
	 */
	void printHeatMap(Tree & tree);

	void sparseBpTable(Tree & tree);

	/**
	 * prints a table with info on each column of the alignment
	 * order: 5' to 3'
	 */
	void columnTable(Tree & tree);

	//consts for realignInterval:
	static string TEMP_FASTA_FILENAME;
	static string TEMP_OUTPUT_FASTA_FILENAME;
	static const string REALIGNER_LOCATION;
	static string JAVA_LOC;
	static const string TCOFFEE_LOC;

	static int minStemLength;

	/**
	 * if printHeaders is true, print header line containing column descriptions
	 * for data generating methods.
	 * Default = true;
	 * TODO: apply to all relevant methods. Currently incomplete.
	 */
	static bool printHeaders;



	//Data structures:
	//TODO: replace gratuitous use of pointers. They should not be necessary...

	vector<string*> seqNames; //vector of sequence names
	vector<string*> alignedSeqs; //vector of aligned sequences (with gaps)
	vector<string*> seqs; //vector of ungapped sequences
	vector<vector<int>* > seq2AlignmentMap; //mapping of sequence positions to alignment positions for each (ungaped) sequence
	vector<vector<int>* > alignment2SeqMap; //mapping of alignment position to sequence position for each (ungapped) sequence
	vector<int> alignedStruct; //vector containing pairing partners of each position in the aligned sequence
							   //positions that are unpaired are set to -1
	vector<vector<int>* > seqStructs; //vector of vectors containing pairing partners of each position for each (ungapped) sequence
									  //positions that are unpaired are set to -1

	vector<int> alignmentLabels; //vector of helix labels for aligned structure
	vector<vector<int>* > helixLabels; //vector of labellings for each (ungapped) sequence
	int helixNumber; //number of true helices

	vector<Helix*> trueHelices; //vector of true helices - does not work with labelHelices2!

	vector<vector<CompetingHelix*>* > competingHelices; //vector of competing helices for each sequence

	vector<vector<vector<StatsWrapper*>* >* > StatsMatrix;
	//Stats matrix:
	//Dimension 1: Sequence
	//Dimension 2: True Helices
	//Dimension 3: Competing Helices

	/**
	 * used for realignInterval... should really be gotten rid of
	 */
	string alignmentName;
protected:

	vector<vector<HelixCore*>* > cores; //deprecated
	vector<InterestingRegion*> interestingRegions; //deprecated

	/**
	 * Log-odds matrices:
	 * look-up table to avoid recalculating log-odds scores for columns
	 * If any entry is > 0, that means it has yet to be calculated
	 */
	vector<vector<double>* > felsDoubles;
	vector<double> felsSingles;

	vector<double> * gapFraction;//stores the fraction of gaps for each column
	vector<double> * seqCons;//stores sequence conservation score for each column

	UTMatrix * consensusBP;

    void readAlignment(string & filename, bool includesStruct = false);
    void readStruct(string& filename);

    /**
     * Sets the structure information to be empty (all unpaired).
     */
    void emptyStruct();

    /*
     * Parses a dot bracket string to fill alignedStruct and seqStructs data structures
     */
    void parseDotBracket(string & dotBracketString);

    /*
     * Populates seqStructs data structure by projecting aligned structure onto individual sequences
     * in the alignment, and making sure the resulting base pairs are valid.
     */
    void fillSeqStructs();

    // bool validBP(char a, char b);
    //validBP moved to Utilities class

    void labelHelices();

    /**
     * Labels each helix with an integer number >0, and prints out a string representation of the
     * labeled sequence.
     */
    void labelHelices2();

    bool isHelixCompeting(int seqIndex, int start_i, int start_j, int length);

    /*
     * adds a helix to a set of helix groups. If it forms a new group, return true.
     * If it is added to an existing group, return false
     */
    bool addHelixToGroups(vector<HelixGroup* > & existingGroups, CompetingHelix & newHelix, int newHelixSeqIndex);

    bool isCompatible(HelixGroup & group, CompetingHelix & helix, int helixSeqIndex);

    void initializeStatsMatrix();

    /**
     * clears all stored vectors containing stats and sequence information
     */
    void clearAll();

	/**
	 * A Clean Competing helix is one in which none of the basepairs
	 * of the helix compete on both sides with a single true helix.
	 */
	bool helixIsClean(CompetingHelix & helix, int seqIndex);

	/*
	 * Calculates the log-likelihood of a competing helix from sequence 'seqIndex', normalized for length
	 */
	double logLikelihood(Tree & tree, CompetingHelix & h, int seqIndex);


	double consensusBPPercent(CompetingHelix & helix, int seqIndex);
	//double noGapCBPPercent(int pos5, int pos3); //not used - need to calculate on a per helix basis
	double noGapCBPPercent(CompetingHelix & helix, int seqIndex);

	/**
	 * returns the fraction of sequences that have a gap at specified alignment position
	 */
	double getGapFraction(unsigned int pos);
	double getGapFraction(CompetingHelix & helix, int seqIndex);


	double getSeqCons(CompetingHelix & helix, int seqIndex);

	double getCovariance(CompetingHelix & helix, int seqIndex);



};

#endif /* ALIGNMENT_H_ */
