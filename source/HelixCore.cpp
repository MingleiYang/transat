/*
 * Copyright: Nicholas P Wiebe (2008-2009) and Irmtraud M Meyer (2008-2009)
 * License: licensed under the GNU General Public License version 3 (GPLv3)
 */

#include "HelixCore.h"
#include <cassert>


HelixCore::HelixCore(int participatingSeqs_, int pos5_, int pos3_):
participatingSeqs(participatingSeqs_), pos5(pos5_), pos3(pos3_){
	length = 1;
}

HelixCore::~HelixCore() {

}

void HelixCore::growOut()
{
	pos5--;
	pos3++;
	length++;
}

void HelixCore::growIn()
{
	length++;
}

string HelixCore::dotBracket(int alignmentLength)
{
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


string HelixCore::printAlignmentAtCore(vector<string*> & alignedSeqs){
	string output = string("");

	for(unsigned int j = 0; j < alignedSeqs.size(); j++){
		for(int i = 0; i < length; i++){
			output += alignedSeqs[j]->at(pos5+i);
		}
		output+= "...";
		for(int i = length-1; i >=0; i--){

			output += alignedSeqs[j]->at(pos3-i);
		}
		output+= "\n";
	}

	for(int i = 0; i < length; i++){
		output += "(";
	}
	output+= "...";
	for(int i = length-1; i >=0; i--){

		output += ")";
	}
	output+= "\n";
	return output;

}




