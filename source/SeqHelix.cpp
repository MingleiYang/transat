/*
 * Copyright: Nicholas P Wiebe (2008-2009) and Irmtraud M Meyer (2008-2009)
 * License: licensed under the GNU General Public License version 3 (GPLv3)
 */

#include "SeqHelix.h"

SeqHelix::SeqHelix(int pos5, int pos3, int length, vector<int>* seq2AlignmentMap_)
: Helix(pos5, pos3, length)
{
	seq2AlignmentMap = seq2AlignmentMap_;

}

SeqHelix::~SeqHelix() {

}

//bool SeqHelix::isSubset(SeqHelix & other)
//{
//
//}
//

