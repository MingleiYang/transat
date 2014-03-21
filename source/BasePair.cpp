/*
 * Copyright: Nicholas P Wiebe (2008-2009) and Irmtraud M Meyer (2008-2009)
 * License: licensed under the GNU General Public License version 3 (GPLv3)
 */

#include "BasePair.h"

BasePair::BasePair(int ahead_, int behind_) {

	ahead = ahead_;
	behind = behind_;

}

BasePair::~BasePair() {

}

bool BasePair::comp(BasePair *a, BasePair *b)
{
	if (a->ahead < b->ahead){
		return true;
	}
	else if (a->ahead > b->ahead){
		return false;
	}
	else{
		return a->behind < b->behind;
	}
}



