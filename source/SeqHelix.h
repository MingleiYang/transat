/*
 * Copyright: Nicholas P Wiebe (2008-2009) and Irmtraud M Meyer (2008-2009)
 * License: licensed under the GNU General Public License version 3 (GPLv3)
 * SeqHelix.h
 * Class for storing Helices that are specific to a certain sequence in the alignment.
 */

#ifndef SEQHELIX_H_
#define SEQHELIX_H_

#include "Helix.h"
#include <vector>

class SeqHelix: public Helix {
public:
	SeqHelix(int pos5, int pos3, int length, vector<int>* seq2AlignmentMap);
	virtual ~SeqHelix();


	bool isSubset(SeqHelix & other);
	bool isEqual(SeqHelix & other);
	bool isSuperSet(SeqHelix & other);

protected:
	/**
	 * seq2AlignmentMap maps sequences positions to positions in the alignment
	 * e.g. seq2AlignmentMap[10] := position in the (gapped) alignment of position 10
	 * in the (ungapped) sequence.
	 * Alignment positions are always greater or equal to sequence positions
	 * positions are indexed starting at 0.
	 */
	vector<int>* seq2AlignmentMap;
};

#endif /* SEQHELIX_H_ */
