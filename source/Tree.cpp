/*
 * Copyright: Nicholas P Wiebe (2008-2009) and Irmtraud M Meyer (2008-2009)
 * License: licensed under the GNU General Public License version 3 (GPLv3)
 */

#include "Tree.h"
#include <iostream>
#include <fstream>
#include <cassert>
#include <sstream>
#include "EvolModel.h"
#include <cmath>
#include "Alignment.h"
#include <cstdlib>
#include "Utilities.h"
#include <limits>

int Tree::leafCount = 1;
bool Tree::nonGapPair = true;
//bool Tree::felsDoubleUnderflow = false;
//bool Tree::felsSingleUnderflow = false;

//setting this variables to true suppresses the underflow warning...
//for distribution, set to true, because it's a confusing error.
bool Tree::felsDoubleUnderflow = true;
bool Tree::felsSingleUnderflow = true;

Tree::Tree(string filename)
{
	ifstream treeFile(filename.c_str());
	string line;
	string newickString = "";

	if(!treeFile.is_open()){
		cerr << "Error: cannot open tree file " << filename << endl;
		exit(-1);
	}

	while(treeFile.good()){
		getline(treeFile, line);
		newickString += line;
	}

	//cout << "Tree input: " << newickString << endl;

	parent = NULL;
	leftChild = NULL;
	rightChild = NULL;

	list<string*> tokens;
	unsigned int split;
	unsigned int start = 0;
	while(start < newickString.length()){
		split = newickString.find_first_of("(),:;", start);
		string * temp;
		if(split > start){
			temp = new string(newickString.substr(start, split-start));
			tokens.push_back(temp);
		}
		temp = new string(newickString.substr(split, 1));
		tokens.push_back(temp);
		start = split+1;
	}



//	cout << "tokenized: ";
//	for(list<string*>::iterator i = tokens.begin(); i != tokens.end(); i++){
//		cout << **i << "|";
//	}
//	cout << endl;

	createChildren(tokens);
	setUpMatrices();

	for(list<string*>::iterator i = tokens.begin(); i != tokens.end(); i++){
		delete (*i);
	}


}

Tree::~Tree() {
	delete leftChild;
	delete rightChild;
}

Tree::Tree(list<string*> newickStringTokens)
{
	createChildren(newickStringTokens);
	setUpMatrices();
}

Tree::Tree(int leaves, double length){
// 1: 1
// 2: 2
// 3: 4
// 4: 6
// 5: 8
// 6: 10
// 7: 12
// 8: 14
// 9:

	parent = NULL;

	if(leaves >= 2){
		double perBranchLength = length / ((leaves-1) * 2);
		int leftLeaves = leaves / 2 + leaves % 2;
		int rightLeaves = leaves / 2;
		leftChild = new Tree(leftLeaves, perBranchLength * (leftLeaves -1) * 2);
		rightChild = new Tree(rightLeaves, perBranchLength * (rightLeaves -1) * 2);

		leftChild->branchLength = perBranchLength;
		rightChild->branchLength = perBranchLength;
		leftChild->parent = this;
		rightChild->parent = this;

		leftChild->setUpMatrices();
		rightChild->setUpMatrices();

		seqName = "";
	}
	else{
		assert(length == 0);
		stringstream ss;
		ss << Tree::leafCount;
		Tree::leafCount++;
		leftChild = NULL;
		rightChild = NULL;
		seqName = ss.str();
	}

}


double Tree::calcFelsDouble(Alignment & a, int pos5, int pos3){
	double likelihood = 0.0;

	//leaf:
	if(leftChild == NULL && rightChild == NULL){
		map<Tree*,unsigned int> leaf2SeqMap = getLeaf2SeqMap(a.seqNames);
		unsigned int seqIndex = leaf2SeqMap.find(this)->second;

		char base5,base3;
		base5 = a.alignedSeqs[seqIndex]->at(pos5);
		base3 = a.alignedSeqs[seqIndex]->at(pos3);


		//convert Ns to gaps:
		if (base5 == 'N' || base5 == 'n'){
			base5 = '-';
		}
		if (base3 == 'N' || base3 == 'n'){
			base3 = '-';
		}

		vector<int> leftInterpreted, rightInterpreted;

		if(nonGapPair && base5 == '-' && base3 != '-'){
			//case: 5' base is a gap, and the other is not
			rightInterpreted =Utilities::interpret(base3);
			leftInterpreted = interpretNonPairingGap(rightInterpreted);
		}
		else if(nonGapPair && base3 == '-' && base5 != '-'){
			//case 3' base is a gap, and the other is not
			leftInterpreted = Utilities::interpret(base5);
			rightInterpreted = interpretNonPairingGap(leftInterpreted);
		}
		else{
			//interpret chars
			leftInterpreted = Utilities::interpret(base5);
			rightInterpreted =Utilities::interpret(base3);

		}

		for(int i = 0; i < 4; i++){
			for(int j = 0; j < 4; j++){
				fels[4 * i + j] = leftInterpreted[i] * rightInterpreted[j];
			}
		}

	}
	//inner node
	else{
		//set fels array for children
		leftChild->calcFelsDouble(a, pos5, pos3);
		rightChild->calcFelsDouble(a, pos5, pos3);

		// then calculate the likelihood based on the likelihood of the two child nodes
		for(int i = 0; i < 16; i++){
			double temp1 = 0.0;
			double temp2 = 0.0;
			for(int j = 0; j < 16; j++){
				temp1 += leftChild->matrix16[i][j] * leftChild->fels[j];
				temp2 += rightChild->matrix16[i][j] * rightChild->fels[j];
			}
			fels[i] = temp1 * temp2;
			if(fels[i] <= 0.0){
				if(!felsDoubleUnderflow){
					cerr << "Warning: calcFelsDouble: underflow problem...\n";
					felsDoubleUnderflow = true;
				}
				fels[i] = numeric_limits<double>::min();
			}
		}
	}
	//root
	if(parent == NULL){
		for(int i = 0; i < 16; i++){
			likelihood += EvolModel::ePiDouble[i] * fels[i];

		}
//assert(likelihood > 0);
	}

	return likelihood;
}

double Tree::calcFelsSingle(Alignment & a, int pos){
	double likelihood = 0.0;

	//leaf
	if(leftChild == NULL && rightChild == NULL){
		map<Tree*,unsigned int> leaf2SeqMap = getLeaf2SeqMap(a.seqNames);
		unsigned int seqIndex = leaf2SeqMap.find(this)->second;

		//interpret char
		vector<int> interpreted = Utilities::interpret(a.alignedSeqs[seqIndex]->at(pos));

		for(int i = 0; i < 4; i++){
			fels[i] = interpreted[i];
			//note: should fels be averaged, so that it sums to 1? will have to think about this...
			//answer = no. See Felsenstein 1981, Extensions section
		}

	}
	//inner node
	else{
		//set fels array for children
		leftChild->calcFelsSingle(a, pos);
		rightChild->calcFelsSingle(a, pos);

		// then calculate the likelihood here
		for(int i = 0; i < 4; i++){
			double temp1 = 0.0;
			double temp2 = 0.0;
			for(int j = 0; j < 4; j++){
				temp1 += leftChild->matrix4[i][j] * leftChild->fels[j];
				temp2 += rightChild->matrix4[i][j] * rightChild->fels[j];
			}
			fels[i] = temp1 * temp2;
			if(fels[i] <= 0.0){
				if(!felsSingleUnderflow){
					cerr << "Warning: calcFelsSingle: underflow problem...\n";
					felsSingleUnderflow = true;
				}
				fels[i] = numeric_limits<double>::min();
			}
		}
	}
	//root
	if(parent == NULL){
		for(int i = 0; i < 4; i++){
			likelihood += EvolModel::ePiSingle[i] * fels[i];
		}

//		assert(likelihood > 0);
	}

	return likelihood;
}

map<Tree*, unsigned int> Tree::getLeaf2SeqMap(vector<string> & SeqNames){
	vector<string*> temp(SeqNames.size());
	for(unsigned int i = 0; i < temp.size(); i++){
		//temp[i] = new string(SeqNames[i]);
		temp[i] = &SeqNames[i];
	}
	map<Tree*, unsigned int> leaf2SeqMap = getLeaf2SeqMap(temp);
//	for(unsigned int i = 0; i < temp.size(); i++){
//		delete temp[i];
//	}
	return leaf2SeqMap;
}

map<Tree*, unsigned int> Tree::getLeaf2SeqMap(vector<string*> & SeqNames)
{
	map<Tree*, unsigned int> out;
	if (leftChild == NULL && rightChild == NULL){
		unsigned int i;
		for(i = 0; i < SeqNames.size(); i++){
			if(seqName.compare(*SeqNames[i]) == 0){
				break;
			}
		}
		if(i >= SeqNames.size()){
			cerr << "Error: tree sequence name not in set of alignment sequences\n";
			cerr << "unrecognized name: " << seqName << endl;
			cerr << "Alignment names:\n";
			for(i = 0; i < SeqNames.size();i++){
				cerr << *SeqNames[i] << endl;
			}
			exit(-1);
		}
		else{
			out.insert(pair<Tree*, unsigned int>(this, i));
		}
	}
	else{
		out = leftChild->getLeaf2SeqMap(SeqNames);
		map<Tree*, unsigned int> rightMap = rightChild->getLeaf2SeqMap(SeqNames);

		out.insert(rightMap.begin(), rightMap.end());
	}

	return out;
}

void Tree::createChildren(list<string*> & tokens)
{
//	cout << "calling createChildren\n";
	//case: inner node
	if ((*tokens.begin())->compare("(") == 0){
//		for(list<string*>::iterator i = tokens.begin(); i != tokens.end(); i++){
//			cout << **i << " $ ";
//		}
//		cout << endl;

		int openParens = 0;
		list<string*>::iterator cursor = ++tokens.begin();


		//find first comma on same level
//		while(cursor != tokens.end() && ((*cursor)->compare(",") != 0 || openParens > 0)){
//
//			if((*cursor)->compare("(") == 0){
//				openParens++;
//			}
//			else if((*cursor)->compare(")") == 0){
//				openParens--;
//			}
//			cursor++;
//		}
//
//		if(cursor == tokens.end()){
//			cerr << "Error: unable to parse tree file.\n";
//			exit(-1);
//		}

		vector<list<string*>::iterator > commas;
		while(cursor != tokens.end()){
			if((*cursor)->compare("(") == 0){
				openParens++;
			}
			else if((*cursor)->compare(")") == 0){
				openParens--;
			}
			else if((*cursor)->compare(",") == 0 && openParens == 0){
				commas.push_back(cursor);
			}
			cursor++;
		}

		assert(!commas.empty());


		cursor = commas[0];

		//cursor is at first comma
		list<string*> leftTokens(++tokens.begin(), cursor);

//		cout << "Left Tokens: ";
//		for(list<string*>::iterator i = leftTokens.begin(); i != leftTokens.end(); i++){
//			cout << **i << "|";
//		}
//		cout << endl;

		leftChild = new Tree(leftTokens);
		leftChild->parent = this;

		cursor++; //skip comma

		list<string*>::iterator rightTokensEnd;
		if ((*tokens.rbegin())->compare(";") == 0){
			branchLength = 0;

			rightTokensEnd = ----tokens.end();

			//some newick trees define a branch length for the root node.
			//(...):0.001; -> with branch length
			//(...); -> without branch length
			if((*rightTokensEnd)->compare(")") != 0){
				//branch length defined... just ignore it
				rightTokensEnd = --------tokens.end();;
				//cout << "-----> " << **rightTokensEnd << endl;
				assert((*rightTokensEnd)->compare(")") == 0);
			}

		}
		else{
			branchLength = atof((*--tokens.end())->c_str());
			rightTokensEnd = ------tokens.end();
		}

		list<string*> rightTokens(cursor, rightTokensEnd);
		if(commas.size() > 1){
			rightTokens.push_front(new string("("));
			rightTokens.push_back(new string(")"));
			rightTokens.push_back(new string(":"));
			rightTokens.push_back(new string("0.0000001"));
			cerr << "Warning: input tree not a binary tree. Binary-izing by adding epsilon-length branches\n";

			//TODO: fix memory leak here... these strings are never deleted.
		}
//		cout << "All Tokens: size=" << tokens.size() << endl;
//		for(list<string*>::iterator i = tokens.begin(); i != tokens.end(); i++){
//			cout << **i << "|";
//		}
//		cout << endl;
//
//		cout << "Right Tokens: size=" << rightTokens.size() << endl;
//		//int count = 0;
//		for(list<string*>::iterator i = rightTokens.begin(); i != rightTokens.end(); i++){
//			//cout << **i << "|"<< count << "|";
//			cout << **i << "|";
//			//count++;
//		}
//		cout << endl;

		rightChild = new Tree(rightTokens);
		rightChild->parent = this;

	}
	//case: leaf
	else{
		leftChild = NULL;
		rightChild = NULL;
		seqName = **tokens.begin();
		branchLength = atof((*--tokens.end())->c_str());
	}
	setUpMatrices();

}



void Tree::setUpMatrices(){
	//code adapted from SimulFold.TreeNode.settingTransitions()

	//initialize matrix4
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){
			matrix4[i][j] = 0.0;
		}
	}

	//set values
	for(int k = 0; k < 4; k++){
		double exponent = exp(EvolModel::eDSingle[k] * branchLength);
		for(int i = 0; i < 4; i++){
			for(int j = 0; j < 4; j++){
				matrix4[i][j] += EvolModel::eVSingle[i][k] * exponent * EvolModel::eWSingle[k][j];
			}
		}
	}

	for(int i = 0; i < 16; i++){
		for(int j = 0; j < 16; j++){
			matrix16[i][j] = 0.0;
		}
	}

	for(int k = 0; k < 16; k++){
		double exponent = exp(EvolModel::eDDouble[k] * branchLength);
		for(int i = 0; i < 16; i++){
			for(int j = 0; j < 16; j++){
				matrix16[i][j] += EvolModel::eVDouble[i][k] * exponent * EvolModel::eWDouble[k][j];
			}
		}
	}

}

string Tree::newickString(){
	stringstream ss;
	//if leaf
	if(leftChild == NULL){
		ss << seqName << ":" << branchLength;

	}
	else{
		//internal node
		ss << "(" << leftChild->newickString() << "," << rightChild->newickString() << ")";

		if(parent != NULL){
			ss << ":" << branchLength;
		}
		else{
			//root - branch length is ignored
			ss << ";";
		}
	}
	//cout << ss.str() << endl;
	return ss.str();
}

vector<int> Tree::interpretNonPairingGap(vector<int> & nongap){
	vector<int> iTable;
	assert(nongap.size() == 4);

	for(int i = 0; i < 4; i++){
		iTable.push_back(1);
	}

	if(nongap[0] == 1){
		//A
		iTable[1] = 0;
	}
	if(nongap[1] == 1){
		//U
		iTable[0] = 0;
		iTable[2] = 0;
	}
	if(nongap[2] == 1){
		//G
		iTable[1] = 0;
		iTable[3] = 0;
	}
	if(nongap[3] == 1){
		//C
		iTable[2] = 0;
	}

	//note: special chars 'b', 'd' and 'k' produce an all-zero iTable... So watch out!
	return iTable;
}


double Tree::totalLength(){
	double length = branchLength;

	if(leftChild != NULL){
		length += leftChild->totalLength();
	}
	if(rightChild != NULL){
		length += rightChild->totalLength();
	}
	return length;
}

vector<string> Tree::getSeqNames(){
	vector<string> names;

	if(leftChild == NULL && rightChild == NULL){
		names.push_back(seqName);
	}
	else{
		vector<string> temp = leftChild->getSeqNames();
		names.insert(names.end(), temp.begin(), temp.end());

		temp = rightChild->getSeqNames();
		names.insert(names.end(), temp.begin(), temp.end());

	}

	return names;
}


