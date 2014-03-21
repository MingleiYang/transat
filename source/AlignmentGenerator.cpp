/*
 * Copyright: Nicholas P Wiebe (2008-2009) and Irmtraud M Meyer (2008-2009)
 * License: licensed under the GNU General Public License version 3 (GPLv3)
 * AlignmentGenerator.cpp
 */

#include "AlignmentGenerator.h"
#include "Utilities.h"
#include "EvolModel.h"
#include <cassert>
#include <iostream>
#include <cstdlib>
#include <stack>
#include <fstream>

const double EPSILON = 0.00001;

AlignmentGenerator::AlignmentGenerator() {


}

AlignmentGenerator::~AlignmentGenerator() {

}

Alignment AlignmentGenerator::makeAlignment(Tree & t, Alignment & a){
	map<Tree*, unsigned int> leaf2SeqMap = t.getLeaf2SeqMap(a.seqNames);
	assert(leaf2SeqMap.size() == a.seqNames.size());

	//initialize new alignment strings
	vector<string> newAlignment;
	vector<string> newAlignmentNames;
	for(unsigned int i = 0; i < a.alignedSeqs.size(); i++){
		newAlignment.push_back("");
		newAlignmentNames.push_back(*a.seqNames[i]);
		for(unsigned int j = 0; j < a.alignedSeqs[i]->length(); j++){
			newAlignment[i] += "X";
		}
	}

	vector<char> columnA, columnB;
	for(unsigned int i = 0; i < a.alignedSeqs.size(); i++){
		columnA.push_back('X');
		columnB.push_back('X');
	}

	double singlePiCDF[4];
	double sum = 0;
	for(int i = 0; i < 4; i++){
		sum += EvolModel::ePiSingle[i];
		singlePiCDF[i] = sum;
	}

	assert(sum >= 1 - EPSILON);
	assert(sum <= 1 + EPSILON);

	singlePiCDF[3] = 1;

	double doublePiCDF[16];
	sum = 0;
	for(int i = 0; i < 16; i++){
		sum += EvolModel::ePiDouble[i];
		doublePiCDF[i] = sum;
	}

	assert(sum >= 1 - EPSILON);
	assert(sum <= 1 + EPSILON);
	doublePiCDF[15] = 1;

	//new seed!
	srand((unsigned)time(0));

	//note: a possible alternative to initializing the
	//root sequence from the prior distribution
	//is to select entries at random from the original alignment
	//TODO: implement the above

	for(unsigned int i = 0; i < a.alignedStruct.size(); i++){
		if (a.alignedStruct[i] == -1){


			double r = (double)rand()/(double)RAND_MAX;
		   	int start_nuc = 0;

		   	while (r > singlePiCDF[start_nuc]){
		   		start_nuc++;
		   	}

		   	fillColumnDF(t.leftChild, leaf2SeqMap, columnA, start_nuc);
		   	fillColumnDF(t.rightChild, leaf2SeqMap, columnA, start_nuc);

		   	for(unsigned int j = 0; j < columnA.size(); j++){
		   		newAlignment[j].at(i) =  columnA[j];
		   	}
		   	//for testing: set columnA back to all Xs, so as to be sure that every
		   	//sequence is being set by fillColumnDF
			for(unsigned int i = 0; i < columnA.size(); i++){
				columnA[i] = 'X';
			}

		}
		else if(a.alignedStruct[i] > (int)i){
			double r = (double)rand()/(double)(RAND_MAX);
		   	int start_nuc = 0;

		   	while (r > doublePiCDF[start_nuc]){
				start_nuc++;
		   	}

		   	fillColumnDF(t.leftChild, leaf2SeqMap, columnA, columnB, start_nuc);
		   	fillColumnDF(t.rightChild, leaf2SeqMap, columnA, columnB, start_nuc);

		   	for(unsigned int j = 0; j < columnA.size(); j++){
		   		newAlignment[j].at(i) =  columnA[j];
		   		newAlignment[j].at(a.alignedStruct[i]) = columnB[j];
			}

		   	//for testing: set columnA and B back to all Xs, so as to be sure that every
		   	//sequence is being set by fillColumnDF
			for(unsigned int i = 0; i < columnA.size(); i++){
				columnA[i] = 'X';
				columnB[i] = 'X';
			}

		}
		//else is 3' base of pair - do nothing
	}

	Alignment out(newAlignmentNames, newAlignment, a.alignedStruct);
	return out;
}

void AlignmentGenerator::fillColumnDF(Tree * node, map<Tree*, unsigned int> & leaf2SeqMap, vector<char> & column, int current){
	double transitionCDF[4];
	double sum = 0;

	for(int i = 0; i < 4; i++){
		sum	 += node->matrix4[current][i];
		transitionCDF[i] = sum;
	}

	assert(sum >= 1 - EPSILON && sum <= 1 + EPSILON);
	transitionCDF[3] = 1;

   	double r = (double)rand()/(double)(RAND_MAX);

   	int new_nuc = 0;
   	while (r > transitionCDF[new_nuc]){
   		new_nuc++;
   	}

   	//if leaf
   	if(node->leftChild == NULL && node->rightChild == NULL){
   		column[leaf2SeqMap.find(node)->second] = Utilities::reverseInterpret(new_nuc);
   	}
   	else{
   		//inner node:
		fillColumnDF(node->leftChild, leaf2SeqMap, column, new_nuc);
		fillColumnDF(node->rightChild, leaf2SeqMap, column, new_nuc);
   	}

}

void AlignmentGenerator::fillColumnDF(Tree * node, map<Tree*, unsigned int> & leaf2SeqMap, vector<char> & columnA, vector<char> & columnB, int current){
	double transitionCDF[16];
	double sum = 0;

	for(int i = 0; i < 16; i++){
		sum	 += node->matrix16[current][i];
		transitionCDF[i] = sum;
	}

	assert(sum >= 1 - EPSILON && sum <= 1 + EPSILON);
	transitionCDF[15] = 1;

   	double r = (double)rand()/(double)(RAND_MAX);

   	int new_nuc = 0;
   	//TODO: implement binary search -> makes this step O(log(N)) rather than O(N) where N is alphabet size
   	while (r > transitionCDF[new_nuc]){
   		new_nuc++;
   	}

   	//if leaf
   	if(node->leftChild == NULL && node->rightChild == NULL){
//column[leaf2SeqMap.find(node)->second] = Utilities::reverseInterpret(new_nuc);
   		int leftNuc = new_nuc / 4;
   		int rightNuc = new_nuc % 4;
   		columnA[leaf2SeqMap.find(node)->second] = Utilities::reverseInterpret(leftNuc);
   		columnB[leaf2SeqMap.find(node)->second] = Utilities::reverseInterpret(rightNuc);
   	}
   	else{
   		//inner node:
		fillColumnDF(node->leftChild, leaf2SeqMap, columnA, columnB, new_nuc);
		fillColumnDF(node->rightChild, leaf2SeqMap, columnA, columnB, new_nuc);
   	}
}

Alignment AlignmentGenerator::makeAlignment(Tree & t, const string & filename){

	//read structure from file
	ifstream structFile(filename.c_str());
	string line;
	string dotBracketStruct = "";
	if(!structFile.is_open()){
		cerr << "Error: cannot open alignment file " << filename << endl;
		exit(-1);
	}
	//skip first line
	getline(structFile, line);
	if(line.at(0) != '>'){
		cerr << "Error: improperly formated file " << filename << "\nExpecting first line to start with '>'\n";
		exit(-1);
	}
	while(structFile.good()){
		getline(structFile, line);
		Utilities::TrimSpaces(line);
		dotBracketStruct += line;
	}

	vector<int> alignedStruct = AlignmentGenerator::parseDotBracket(dotBracketStruct);


	vector<string> newAlignmentNames = t.getSeqNames();

	int alignmentSize = newAlignmentNames.size();
	int alignmentLength = alignedStruct.size();

	map<Tree*, unsigned int> leaf2SeqMap = t.getLeaf2SeqMap(newAlignmentNames);
	assert(leaf2SeqMap.size() == newAlignmentNames.size());

	vector<string> newAlignment(alignmentSize, string(alignmentLength, 'X'));


	vector<char> columnA(alignmentSize, 'X');
	vector<char> columnB(alignmentSize, 'X');

	double singlePiCDF[4];
	double sum = 0;
	for(int i = 0; i < 4; i++){
		sum += EvolModel::ePiSingle[i];
		singlePiCDF[i] = sum;
	}

	assert(sum >= 1 - EPSILON);
	assert(sum <= 1 + EPSILON);

	singlePiCDF[3] = 1;

	double doublePiCDF[16];
	sum = 0;
	for(int i = 0; i < 16; i++){
		sum += EvolModel::ePiDouble[i];
		doublePiCDF[i] = sum;
	}

	assert(sum >= 1 - EPSILON);
	assert(sum <= 1 + EPSILON);
	doublePiCDF[15] = 1;

	//new seed!
	srand((unsigned)time(0));

	//note: a possible alternative to initializing the
	//root sequence from the prior distribution
	//is to select entries at random from the original alignment
	//TODO: implement the above

	for(int i = 0; i < alignmentLength; i++){
		if (alignedStruct[i] == -1){


			double r = (double)rand()/(double)RAND_MAX;
			int start_nuc = 0;

			while (r > singlePiCDF[start_nuc]){
				start_nuc++;
			}

			fillColumnDF(t.leftChild, leaf2SeqMap, columnA, start_nuc);
			fillColumnDF(t.rightChild, leaf2SeqMap, columnA, start_nuc);

			for(unsigned int j = 0; j < columnA.size(); j++){
				newAlignment[j].at(i) =  columnA[j];
			}
			//for testing: set columnA back to all Xs, so as to be sure that every
			//sequence is being set by fillColumnDF
			for(unsigned int i = 0; i < columnA.size(); i++){
				columnA[i] = 'X';
			}

		}
		else if(alignedStruct[i] > (int)i){
			double r = (double)rand()/(double)(RAND_MAX);
			int start_nuc = 0;

			while (r > doublePiCDF[start_nuc]){
				start_nuc++;
			}

			fillColumnDF(t.leftChild, leaf2SeqMap, columnA, columnB, start_nuc);
			fillColumnDF(t.rightChild, leaf2SeqMap, columnA, columnB, start_nuc);

			for(unsigned int j = 0; j < columnA.size(); j++){
				newAlignment[j].at(i) =  columnA[j];
				newAlignment[j].at(alignedStruct[i]) = columnB[j];
			}

			//for testing: set columnA and B back to all Xs, so as to be sure that every
			//sequence is being set by fillColumnDF
			for(unsigned int i = 0; i < columnA.size(); i++){
				columnA[i] = 'X';
				columnB[i] = 'X';
			}

		}
		//else is 3' base of pair - do nothing
	}

	Alignment out(newAlignmentNames, newAlignment, alignedStruct);
	return out;
}


vector<int> AlignmentGenerator::parseDotBracket(string & dotBracketString)
{
	//TODO: remove duplicate function in Alignment class
	//cout << dotBracketString.length() << endl;
	stack<int> parenStack, bracketStack, ltgtStack, braceStack, aStack, bStack, cStack, dStack;;

	vector<int> alignedStruct(dotBracketString.length(), -1);

	int pairingPos;

	for(int i = 0; i < (int)dotBracketString.length(); i++){
		char a = dotBracketString.at(i);
		pairingPos = -1;
		switch(a){
		case '(':
			parenStack.push(i);
			break;
		case '<':
			ltgtStack.push(i);
			break;
		case '[':
			bracketStack.push(i);
			break;
		case '{':
			braceStack.push(i);
			break;
		case 'A':
			aStack.push(i);
			break;
		case 'B':
			bStack.push(i);
			break;
		case 'C':
			cStack.push(i);
			break;
		case 'D':
			dStack.push(i);
			break;
		case ')':
			assert(!parenStack.empty());
			pairingPos = parenStack.top();
			parenStack.pop();
			break;
		case '>':
			assert(!ltgtStack.empty());
			pairingPos = ltgtStack.top();
			ltgtStack.pop();
			break;
		case ']':
			assert(!bracketStack.empty());
			pairingPos = bracketStack.top();
			bracketStack.pop();
			break;
		case '}':
			assert(!braceStack.empty());
			pairingPos = braceStack.top();
			braceStack.pop();
			break;
		case 'a':
			assert(!aStack.empty());
			pairingPos = aStack.top();
			aStack.pop();
			break;
		case 'b':
			assert(!bStack.empty());
			pairingPos = bStack.top();
			bStack.pop();
			break;
		case 'c':
			assert(!cStack.empty());
			pairingPos = cStack.top();
			cStack.pop();
			break;
		case 'd':
			assert(!dStack.empty());
			pairingPos = dStack.top();
			dStack.pop();
			break;
		case '.':
			//do nothing - pairing pos = -1 --> not paired
			break;
		default:
			cerr << "Error: unrecognized character '" << a << "' in dot bracket string\n";
			exit(-1);

		}

		alignedStruct[i] = pairingPos;
		if(pairingPos >= 0){
			alignedStruct[pairingPos] = i;
		}

	}

	assert(parenStack.empty());
	assert(bracketStack.empty());
	assert(ltgtStack.empty());
	assert(braceStack.empty());
	assert(aStack.empty());
	assert(bStack.empty());
	assert(cStack.empty());
	assert(dStack.empty());

	return alignedStruct;

}

