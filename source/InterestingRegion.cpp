/*
 * Copyright: Nicholas P Wiebe (2008-2009) and Irmtraud M Meyer (2008-2009)
 * License: licensed under the GNU General Public License version 3 (GPLv3)
 */

#include "InterestingRegion.h"
#include <iostream>


InterestingRegion::InterestingRegion(int start_, int end_, int competedHelix_, double structCons_, double seqCons_) {

	start = start_;
	end = end_;
	competedHelix = competedHelix_;
	seqConservation = seqCons_;
	structConservation = structCons_;
}

InterestingRegion::~InterestingRegion() {

}

double InterestingRegion::score()
{
	return structConservation / seqConservation;
}

void InterestingRegion::print()
{
	cout << "score:\t" << score() << "\tseqCons:\t" << seqConservation << "\tstructCons:\t" << structConservation << "\tstart:\t" << start << "\tend:\t" << end << "\thelix:\t" << competedHelix << endl;
}

bool InterestingRegion::comp(InterestingRegion *a, InterestingRegion *b)
{
	return a->score() < b->score();
}




