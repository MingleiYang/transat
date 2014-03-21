/*
 * Copyright: Nicholas P Wiebe (2008-2009) and Irmtraud M Meyer (2008-2009)
 * License: licensed under the GNU General Public License version 3 (GPLv3)
 * CompetingHelix.h
 * Stores information about a single competing helix.
 * Positions are with respect to the sequence in which the competing helix is found
 * (i.e. they are not alignment positions)
 */

#ifndef COMPETINGHELIX_H_
#define COMPETINGHELIX_H_

#include <string>
#include <vector>
#include "StatsWrapper.h"

using namespace std;

class CompetingHelix {
public:
	//deprecated:
	CompetingHelix(int _start5, int _end5, int _start3, int _end3, StatsWrapper _stat);
	//use this one:
	CompetingHelix(int _start5, int _end5, int _start3, int _end3);
	virtual ~CompetingHelix();


	//helix looks like this:
	//------------((((((((((---------------------))))))))))--------
	//    5'start ^  5'end ^              3'start^   3'end^
	int start5, end5, start3, end3;


	StatsWrapper stat; //don't use this -> use StatsMatrix in Alignment class instead

	string printHelix(int iL);
	string printHelix2(int iL);
	string printRevHelix(int iL);

	string printAlignedHelix(vector<int> &seq2AlignmentMap, int alignmentLength);


	unsigned int isPosCompeting(int iPos);

	double midpoint();

	double freeEnergy(const double eStack[4][4][4][4], string & seq);



};

#endif /* COMPETINGHELIX_H_ */
