/*
 * Copyright: Nicholas P Wiebe (2008-2009) and Irmtraud M Meyer (2008-2009)
 * License: licensed under the GNU General Public License version 3 (GPLv3)
 * AlignedHelix.h
 * Class for storing a helix as a set of base-pairs
 */

#ifndef ALIGNEDHELIX_H_
#define ALIGNEDHELIX_H_

#include <vector>
#include <set>
#include <string>
#include "Helix.h"
#include "Tree.h"

using namespace std;

class Alignment;

class AlignedHelix {
public:
	AlignedHelix(const Helix & helix, int seqIndex, vector<int> & seq2AlignmentMap);
	AlignedHelix(int outerBpPos5, int outerBpPos3, int length);
	virtual ~AlignedHelix();

	string dotBracket(unsigned int alignmentLength);

	/**
	 * methods to collect stats on a give helix
	 */
	double pairedLikelihood(Alignment & a, Tree & t);
	double unpairdLikelihood(Alignment & a, Tree & t);
	double logLikeRatio(Alignment & a, Tree & t);
	double canonicalBP(Alignment & a);
	double covariance(Alignment & a);
	double conservation(Alignment & a);


	/**
	 * tries to insert the given helix
	 * If the helix aligns to exactly the same set of bps as this helix, add this sequence
	 * to the list of sequences in which this alignedHelix appears and return true.
	 * If not, do nothing and return false.
	 */
	bool insert(const Helix & helix, int seqIndex, vector<int> & seq2AlignmentMap);

	/**
	 * returns true if helix is found in alignment consensus structure, else false
	 *
	 * exact definition: all bps of helix are in consensus structure, and consensus
	 * structure does not have bps inner to or outer to any of the bps in the helix
	 * that are not also in the helix.
	 *
	 * example:
	 * ...((((...))))... consensus
	 * ...((((...))))... helix
	 * returns true
	 *
	 * ...((((...))))... consensus
	 * ..(((((...))))).. helix
	 * returns false
	 *
	 * ...((((...))))... consensus
	 * ....(((...))).... helix
	 * returns false
	 *
	 * ...((((...))))... consensus
	 * ...((.(...).))... helix
	 * returns false
	 */
	bool isConsensusHelix(Alignment & a);

	/**
	 * returns true if any of the bps of the helix are in the consensus structure
	 */
	bool isPartialConsensusHelix(Alignment & a);

	/**
	 * returns the number of bps this helix shares with the consensus structure
	 */
	int consensusBps(Alignment & a);

	/**
	 * returns true if the helix is a competing helix in any of the sequences in which it
	 * appears.
	 */
	bool isCompetingHelix(Alignment & a);

	/**
	 * returns the average midpoint of bps in the helix (with respect to the alignment as a whole)
	 */
	double midpoint();

	/**
	 * returns the average midpoint of bps in the helix (with respect to the alignment as a whole),
	 * normalized for the alignment length
	 */
	double midpoint(Alignment & a);

	/**
	 * returns the number of bps in the helix
	 */
	int length();

	/**
	 * returns the number of sequences in which this exact helix appears
	 * note: this does not count sub- or supersets of the helix
	 */
	int appearances();

	/**
	 * returns a string with the positions of the bps, separated by commas
	 * example:
	 * 12:32,11:31,10:29
	 */
	string bpsString();

	pair<double,double> cisTransScore(Alignment & a);

	void competeScore(Alignment & a, double & cis5, double & cis3, double & trans5, double & trans3, double & mid5, double & mid3);
	//double trans(Alignment & a);

	/**
	 * vector containing pairs of integers representing the position of each bp of the helix.
	 * First in pair is 5' position
	 * Second in pair is 3' position
	 * Positions are with respect to the alignment.
	 * Bps in this vector are ordered outer to inner.
	 *
	 */
	vector<pair<int, int> > bps; //could be a linked list instead... not sure it would make a difference


private:
	//TODO: this class should probably store a pointer to
	//the alignment in which it exists in, so that the methods
	//don't need to take as parameters the alignment


	/**
	 * set containing the alignment indexes of the sequences in which this helix appears
	 */
	set<int> appearsIn;

	/**
	 * list of consensus helices that this aligned helix competes with
	 */
	list<AlignedHelix*> competesWith;



};

#endif /* ALIGNEDHELIX_H_ */
