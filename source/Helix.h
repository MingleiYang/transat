/*
 * Copyright: Nicholas P Wiebe (2008-2009) and Irmtraud M Meyer (2008-2009)
 * License: licensed under the GNU General Public License version 3 (GPLv3)
 * Helix.h
 *
 * Class for storing helices of various kinds
 * Stores start point, end point, and length;
 *
 * TODO make CompetingHelix a subclass of this class
 */

#ifndef HELIX_H_
#define HELIX_H_

#include <string>
#include <vector>

using namespace std;

class Helix {
public:

	Helix(int pos5 = 0, int pos3 = 0, int length = -1);
	virtual ~Helix();


	int pos5; //5'-most base
	int pos3; //pairing partner of 5'-most base
	int length;
	// example:
	// --------(((---)))----
	// ----pos5^---pos3^----
	//length = 3

	/*
	 * calculate the midpoint of the helix, divided by the normalization factor
	 * To get midpoint normalized for sequence length, set the normalization factor
	 * to the sequence length.
	 */
	double midpoint(double normalize = 1.0);

	double freeEnergy(const double eStack[4][4][4][4], string & seq);

	string printHelix(int AlignmentLength);


};

#endif /* HELIX_H_ */
