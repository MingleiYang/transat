/*
 * Copyright: Nicholas P Wiebe (2008-2009) and Irmtraud M Meyer (2008-2009)
 * License: licensed under the GNU General Public License version 3 (GPLv3)
 * AlignmentGenerator.h
 */

#ifndef ALIGNMENTGENERATOR_H_
#define ALIGNMENTGENERATOR_H_

#include "Tree.h"
#include "Alignment.h"
#include <vector>
#include <string>


using namespace std;

class AlignmentGenerator {
public:
	AlignmentGenerator();
	virtual ~AlignmentGenerator();

	/**
	 * Generates a new alignment from a given tree t and the consensus structure from
	 * alignment a, using reverse felsenstein algorithm.
	 *
	 * these should probably be made constructors, so that this is more of a true factory
	 */
	Alignment makeAlignment(Tree & t, Alignment & a);
	Alignment makeAlignment(Tree & t, const string & filename);

private:

	/**
	 * recursive method for Depth-first traversal of tree, filling the column vector with randomly generated
	 * non-paired alignment column
	 *
	 * @param current: nuceotide at the start of the branch defined by the current node
	 */
	void fillColumnDF(Tree * node, map<Tree*, unsigned int> & leaf2SeqMap, vector<char> & column, int current);
	void fillColumnDF(Tree * node, map<Tree*, unsigned int> & leaf2SeqMap, vector<char> & columnA, vector<char> & columnB, int current);

	static vector<int> parseDotBracket(string & dotBracketString);
};



#endif /* ALIGNMENTGENERATOR_H_ */
