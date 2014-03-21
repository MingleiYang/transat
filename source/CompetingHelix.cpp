/*
 * Copyright: Nicholas P Wiebe (2008-2009) and Irmtraud M Meyer (2008-2009)
 * License: licensed under the GNU General Public License version 3 (GPLv3)
 */

#include "CompetingHelix.h"
#include <iostream>
#include <cassert>
#include <sstream>
#include "Utilities.h"

CompetingHelix::CompetingHelix(int _start5, int _end5, int _start3, int _end3, StatsWrapper _stat)
: start5(_start5), end5(_end5), start3(_start3), end3(_end3), stat(_stat) {
	assert(start5 <= end5);
	assert(end5 < start3);
	assert(start3 <= end3);
}
/*
CompetingHelix::CompetingHelix(int _start5, int _end5, int _start3, int _end3, StatsWrapper _stat)
: start5(_start5), end5(_end5), start3(_start3), end3(_end3), statsWrapper(_stat) {
	assert(start5 <= end5);
	assert(end5 < start3);
	assert(start3 <= end3);
} */
CompetingHelix::CompetingHelix(int _start5, int _end5, int _start3, int _end3)
: start5(_start5), end5(_end5), start3(_start3), end3(_end3) {
	assert(start5 <= end5);
	assert(end5 < start3);
	assert(start3 <= end3);
}

CompetingHelix::~CompetingHelix() {
//	if (stat != NULL) delete stat;
}

string type_to_string (HelixType type) {
	string str;

	switch (type) {
	case NONE:
		str = "none";
		break;
	case CIS5:
		str = "5'cis";
		break;
	case TRANS5:
		str = "5'trans";
		break;
	case CIS3:
		str = "3'cis";
		break;
	case TRANS3:
		str = "3'trans";
		break;
	case IN:
		str = "in";
		break;
	default:
		str = "none";
	}

	return str;
}

//returns a string in dot bracket notation representing the competing helix
string CompetingHelix::printHelix2(int iL){

	string output = "";

	int i = 0;
	for (i = 0; i < start5; i++){
		output += ".";
	}

	for (;i <= end5; i++){
		output += "(";
	}

	for (;i < start3; i++){
		output +=".";
	}
	for (;i <= end3; i++){
		output += ")";
	}

	for (;i < iL; i++){
		output += ".";
	}

	return output;
}

string CompetingHelix::printAlignedHelix(vector<int> & seq2AlignmentMap, int alignmentLength)
{
	//int alignedStart5, alignedEnd5, alignedStart3, alignedEnd3;

	int alignedStart5 = seq2AlignmentMap[start5];
//	alignedEnd5 = seq2AlignmentMap[end5];
//	alignedStart3 = seq2AlignmentMap[start3];
	int alignedEnd3 = seq2AlignmentMap[end3];

	string output = "";
	for(int j = 0; j < alignmentLength; j++){
		output += ".";
	}

	for(int j = 0; j < end5-start5+1; j++){
		alignedStart5 = seq2AlignmentMap[start5 + j];
		alignedEnd3 = seq2AlignmentMap[end3 - j];
		output[alignedStart5] = '(';
		output[alignedEnd3] = ')';
	}

//	int i = 0;
//	for (i = 0; i < alignedStart5; i++){
//		output += ".";
//	}
//
//	for (;i <= alignedEnd5; i++){
//		output += "(";
//	}
//
//	for (;i < alignedStart3; i++){
//		output +=".";
//	}
//	for (;i <= alignedEnd3; i++){
//		output += ")";
//	}
//
//	for (;i < alignmentLength; i++){
//		output += ".";
//	}

	return output;

}

// returns 1 if the position is in a competing helix
unsigned int CompetingHelix::isPosCompeting(int iPos){

	if( ((iPos >= start5) && (iPos <= end5))
		||
		((iPos >= start3) && (iPos <= end3)) ) {
			return 1;
	} else {
			return 0;
	}
}



string CompetingHelix::printHelix(int iL)
{
	assert(end3 < iL);

	string output = string("");

	// print out the type of helix and its stat value at the beginning
	// of the line for this competing helix
	stringstream out;
	out << type_to_string (stat.type) << " pval: " << stat.pStat << endl;
	output += out.str();
	//output += " pval: ";
	//output += dtoa(stat.pStat);
	//output += " ";

	int i = 0;
	for (i = 0; i < start5; i++){
		output += ".";
	}

	for (;i <= end5; i++){
		output += "(";
	}

	for (;i < start3; i++){
		output +=".";
	}
	for (;i <= end3; i++){
		output += ")";
	}

	for (;i < iL; i++){
		output += ".";
	}

	return output;
}

string CompetingHelix::printRevHelix(int iL){

	assert(end3 < iL);

	string output = string("");

	int i = 0;
	for (i = 0; i < start5; i++){
		output = "." + output;
	}

	for (;i <= end5; i++){
		output = ")" + output;
	}

	for (;i < start3; i++){
		output = "." + output;
	}
	for (;i <= end3; i++){
		output = "(" + output;
	}

	for (;i < iL; i++){
		output = "." + output;
	}

	stringstream out;
	out << type_to_string (stat.type) << " pval: " << stat.pStat << endl;
	output = out.str() + output;

	return output;
}

double CompetingHelix::midpoint()
{
	return (end3+start5) / 2.0;
}

double CompetingHelix::freeEnergy(const double eStack[4][4][4][4], string & seq)
{
	double eE = 0;
	//loop over bp 'stacks' (adjacent bps)
	//k=0:length-2
	for (int k = 0; k < end5 - start5; k++) {

		int i_2 = start5 + k;
		int j_2 = end3 - k;

		eE += eStack[Utilities::iAt(seq[i_2])][Utilities::iAt(seq[i_2 + 1])][Utilities::iAt(seq[j_2 - 1])][Utilities::iAt(seq[j_2])];
	}
	return eE;
}




