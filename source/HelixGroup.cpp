/*
 * Copyright: Nicholas P Wiebe (2008-2009) and Irmtraud M Meyer (2008-2009)
 * License: licensed under the GNU General Public License version 3 (GPLv3)
 */

#include "HelixGroup.h"
#include <cassert>


HelixGroup::HelixGroup(int numberOfSequences, int alignmentLength, int competedHelix_)
{
	competedHelix = competedHelix_;
	for(int i = 0; i < numberOfSequences; i++){
		group.push_back(NULL);
	}

	//for(int i = 0; i < alignmentLength; i++){
	//	pairingPartners.push_back(-1);
	//}
	groupSize = 0;

}

HelixGroup::~HelixGroup() {
	for(unsigned int i = 0; i < group.size(); i++){
		if(group.at(i) != NULL){
			delete group.at(i);
		}
	}
	group.clear();
}

void HelixGroup::addCompetingHelix(int seqIndex, CompetingHelix* helix, vector<vector<int>* > & seq2AlignmentMap)
{
	assert(seqIndex < (int) group.size());

	if(group[seqIndex] == NULL){
		groupSize++;
	}

	group[seqIndex] = helix;

	int start5 = helix->start5;
	int end5 = helix->end5;
	//int start3 = helix->start3;
	int end3 = helix->end3;

	if (groupSize > 1){
		while(start5 <= end5){
			//int start5align = seq2AlignmentMap[seqIndex]->at(start5);
			//int end3align = seq2AlignmentMap[seqIndex]->at(end3);


			start5++;
			end3--;
		}
	}
	else{
		//first helix in group

		while(start5 <= end5){
//			int start5align = seq2AlignmentMap[seqIndex]->at(start5);
//			int end3align = seq2AlignmentMap[seqIndex]->at(end3);

//			corePairs.insert(start5align, end3align);

			start5++;
			end3--;
		}
	}
}


