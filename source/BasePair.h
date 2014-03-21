/*
 * Copyright: Nicholas P Wiebe (2008-2009) and Irmtraud M Meyer (2008-2009)
 * License: licensed under the GNU General Public License version 3 (GPLv3)
 * BasePair.h
 */

#ifndef BASEPAIR_H_
#define BASEPAIR_H_

class BasePair {
public:
	BasePair(int ahead, int behind);
	virtual ~BasePair();

	int ahead;
	int behind;

	static bool comp(BasePair* a, BasePair* b);
};

#endif /* BASEPAIR_H_ */
