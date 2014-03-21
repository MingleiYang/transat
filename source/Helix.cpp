/*
 * Copyright: Nicholas P Wiebe (2008-2009) and Irmtraud M Meyer (2008-2009)
 * License: licensed under the GNU General Public License version 3 (GPLv3)
 */

#include "Helix.h"
#include "Utilities.h"
#include <cassert>

Helix::Helix(int pos5_, int pos3_, int length_) {
	pos5 = pos5_;
	pos3 = pos3_;
	length = length_;
}

Helix::~Helix() {

}

double Helix::midpoint(double normalize)
{
	return (pos3+pos5) / (normalize * 2.0);
}

double Helix::freeEnergy(const double eStack[4][4][4][4], string & seq)
{
	double eE = 0;
	//loop over bp 'stacks' (adjacent bps)
	//k=0:length-2
	for (int k = 0; k < length-1; k++) {

		int i_2 = pos5 + k;
		int j_2 = pos3 - k;
		if(seq[i_2] != '-' && seq[i_2 + 1] != '-' && seq[j_2] != '-' && seq[j_2-1] != '-'){
			eE += eStack[Utilities::iAt(seq[i_2])][Utilities::iAt(seq[i_2 + 1])][Utilities::iAt(seq[j_2 - 1])][Utilities::iAt(seq[j_2])];
		}
	}
	return eE;
}

string Helix::printHelix(int alignmentLength){
	assert(pos3 < alignmentLength);

	string output = string("");
	int i;
	for(i = 0; i < pos5; i++){
		output += ".";
	}
	for(; i < pos5 + length; i++){
		output += "(";
	}
	for(; i < pos3 - length +1; i++){
		output += ".";
	}
	for(; i < pos3 + 1; i++){
		output += ")";
	}
	for(; i < alignmentLength; i++){
		output += ".";
	}
	assert((int)output.length() == alignmentLength);

	return output;

}

