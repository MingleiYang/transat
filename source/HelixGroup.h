/*
 * Copyright: Nicholas P Wiebe (2008-2009) and Irmtraud M Meyer (2008-2009)
 * License: licensed under the GNU General Public License version 3 (GPLv3)
 * HelixGroup.h
 */

#ifndef HELIXGROUP_H_
#define HELIXGROUP_H_

#include <vector>
#include <map>
#include "CompetingHelix.h"

using namespace std;

class HelixGroup {
public:
	HelixGroup(int numberOfSequences, int alignmentLength, int competedHelix);
	virtual ~HelixGroup();

	void addCompetingHelix(int helixIndex, CompetingHelix* helix, vector<vector<int>* > & seq2AlignmentMap);


	//getters and setters:
	int getCompetedHelix() const
    {
        return competedHelix;
    }

    vector<CompetingHelix*> getGroup() const
    {
        return group;
    }

    int getGroupSize() const
    {
        return groupSize;
    }

private:
	int competedHelix; //index of helix competed with
	vector<CompetingHelix* > group; //this vector contains one entry for each sequence
									  //in the alignment. Sequences without a competing helix
									  //in the group are assigned null
	vector<int> pairingPartners;
	map<int, int> corePairs;

	int groupSize;


};

#endif /* HELIXGROUP_H_ */
