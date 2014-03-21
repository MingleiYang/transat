/*
 * Copyright: Nicholas P Wiebe (2008-2009) and Irmtraud M Meyer (2008-2009)
 * License: licensed under the GNU General Public License version 3 (GPLv3)
 * StatsWrapper.h
 * Wrapper for iStats, eStats arrays used to calculate cis and trans values.
 */

#ifndef STATSWRAPPER_H_
#define STATSWRAPPER_H_

enum HelixType { NONE, CIS5, TRANS5, CIS3, TRANS3, IN }; //deprecated
//note: a competing helix may compete in more that one way, so this is not necessarily
//a good solution. For instance, if the true structure is:
// ....(((...))).(((...)))...
//with competing helix:
// ...........(((((.....)))))
//this helix competes in two different ways with two true helices. -njwiebe

class StatsWrapper {
public:
	StatsWrapper();
	StatsWrapper(const StatsWrapper &other);
	virtual ~StatsWrapper();
	int* iStats;
	double* eStats;
	int* revIStats;
	double* revEStats;

	double pStat;
	double gStat;
	double revgStat;
	HelixType type;
	HelixType revtype;

	double cis5();
	double trans5();
	double cis3();
	double trans3();
	double mid5(); // (i^, c, i) 5' mid --> competing base is between i and i^, and 5' of the base it pairs with in the competing helix
	double mid3(); // (i, c, i^) 3' mid --> competing base is between i and i^, and 3' of the base it pairs with in the competing helix
	double trans();
	double cis();
	double mid();

	bool isZero();

	//operator overloading, for easy summing of several StatsWrapper objects
	//based on examples in http://www.cs.caltech.edu/courses/cs11/material/cpp/donnie/cpp-ops.html
	StatsWrapper & operator+=(const StatsWrapper &rhs);
	const StatsWrapper operator+(const StatsWrapper &other) const;
};

#endif /* STATSWRAPPER_H_ */
