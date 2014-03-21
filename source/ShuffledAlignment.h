/*
 * Copyright: Nicholas P Wiebe (2008-2009) and Irmtraud M Meyer (2008-2009)
 * License: licensed under the GNU General Public License version 3 (GPLv3)
 * shuffledAlignment.h
 * Creates a randomized alignment from an existing alignment
 */

#ifndef SHUFFLEDALIGNMENT_H_
#define SHUFFLEDALIGNMENT_H_

#include "Alignment.h"
#include "Tree.h"
#include <string>

using namespace std;

class ShuffledAlignment: public Alignment {
public:

	/**
	 * Creates a shuffled alignment from the given alignment a
	 * Original structure is discarded. Shuffling is performed with
	 * RNAz shuffling perl script.
	 */
	ShuffledAlignment(const Alignment & a);
	/*
	 * Create shuffled alignment from the given alignment a
	 * Shuffle alignment, keeping only columns in true helix 'th' fixed
	 */
	ShuffledAlignment(const Alignment & a, int th, string treeFile = "");

	virtual ~ShuffledAlignment();

//static const string RNAZ_SHUFFLER_LOC;
protected:

	/**
	 * shuffles entire alignment using RNAz shuffler (rnazRandomizeAln.pl)
	 * Uses level 1 shuffling (columns binned by mean pairwise identity,
	 * rounded to the nearest 10%).
	 */
	void shuffleRNAz();

	/**
	 * shuffles specified columns using RNAz shuffler (rnazRandomizeAln.pl)
	 * Uses level 1 shuffling (columns binned by mean pairwise identity,
	 * rounded to the nearest 10%).
	 */
	void shuffleRNAz(const vector<int> & columns);

	/**
	 * reads a file in clustalW format into alignedSeqs
	 * TODO: move to Alignment Class (no reason to be here and not there).
	 */
	void parseClustalWFile(const string & filename);

	/**
	 * reads a file in clustalW format into alignedSeqs, slotting alignment into
	 * specified columns
	 */
	void parseClustalWFile(const string & filename, const vector<int> & columns);


};



#endif /* SHUFFLEDALIGNMENT_H_ */
