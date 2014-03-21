/*
 * Copyright: Nicholas P Wiebe (2008-2009) and Irmtraud M Meyer (2008-2009)
 * License: licensed under the GNU General Public License version 3 (GPLv3)
 * InterestingRegion.h
 */

#ifndef INTERESTINGREGION_H_
#define INTERESTINGREGION_H_

using namespace std;

class InterestingRegion {
public:
	InterestingRegion(int start_, int end_, int competedHelix_, double structCons_, double seqCons_);
	virtual ~InterestingRegion();

	double score(); // = seqConservation / structConservation
	double seqConservation;
	double structConservation;

	int competedHelix;
	int start, end;

	void print();

	static bool comp(InterestingRegion* a, InterestingRegion* b);
};

#endif /* INTERESTINGREGION_H_ */
