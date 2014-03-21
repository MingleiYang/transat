/*
 * Copyright: Nicholas P Wiebe (2008-2009) and Irmtraud M Meyer (2008-2009)
 * License: licensed under the GNU General Public License version 3 (GPLv3)
 */

#include "Utilities.h"
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <string>

Utilities::Utilities() {

}

int Utilities::iAt(char cC)
{

	if (cC == 'A') {
		return (0);
	} else if (cC == 'C') {
		return (1);
	} else if (cC == 'G') {
		return (2);
	} else if (cC == 'U') {
		return (3);
	} else {
		cerr << "Hey Huston, we have a problem... " << cC << endl;
		assert(false);
	}

	return (-100);

}

bool Utilities::iCheckLongHelix(vector<int> & Bps, int pos)
{

	bool iJ = false;

	if ((pos + 2 < (int)Bps.size()) && (Bps[pos + 1] == Bps[pos] - 1) && (Bps[pos + 2] == Bps[pos] - 2)) {
		iJ = true;
	}
	else if ((pos >= 2) && (Bps[pos - 1] == Bps[pos] + 1) && (Bps[pos - 2] == Bps[pos] + 2)) {
		iJ = true;
	}
	else if ((pos >= 1) && (pos + 1 < (int)Bps.size()) && (Bps[pos + 1] == Bps[pos] - 1) && (Bps[pos - 1] == Bps[pos] + 1)) {
		iJ = true;
	}

	return iJ;
}

int Utilities::iCheckLongHelix(const int* iPosition, const int iI)
{
	int iJ;
	iJ = 0;
	if ((iPosition[iI + 1] == iPosition[iI] - 1) && (iPosition[iI + 2]
			== iPosition[iI] - 2)) {
		iJ = 1;
	}
	if ((iPosition[iI - 1] == iPosition[iI] + 1) && (iPosition[iI - 2]
			== iPosition[iI] + 2)) {
		iJ = 1;
	}
	if ((iPosition[iI + 1] == iPosition[iI] - 1) && (iPosition[iI - 1]
			== iPosition[iI] + 1)) {
		iJ = 1;
	}

	return (iJ);
}

void Utilities::ReadStack(double eStack[4][4][4][4]){
	FILE *fStack;
	int i, j, k, l;

	fStack = fopen("stack_with_zeros.37", "r");
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			for (l = 0; l < 4; l++) {
				for (k = 0; k < 4; k++) {
					fscanf(fStack, "%lf", &eStack[i][j][k][l]);
					//cout << eStack[i][j][k][l] << endl;
				}
			}
		}
	}
	fclose(fStack);
}

int Utilities::lowerBound(vector<double> & v, int begin, int end, double value){
	int diff = begin-end;
	if (diff == 0){
		return begin;
	}
	else{
		int mid = (begin + end) / 2;
		assert(mid < (int)v.size());
		if(v[mid] < value){
			return Utilities::lowerBound(v, mid+1, end, value);
		}
		else{
			return Utilities::lowerBound(v, begin, mid, value);
		}
	}
}

Utilities::~Utilities() {
}

vector<int> Utilities::interpret(char cC)
{
	vector<int> iTable;
	for(int i = 0; i < 4; i++){
		iTable.push_back(0);
	}

	if(cC == 'a' || cC == 'A'){
		iTable[0] = 1;
		iTable[1] = 0;
		iTable[2] = 0;
		iTable[3] = 0;
	}
	else if(cC == 'u' || cC == 'U' || cC == 't' || cC == 'T'){
		iTable[0] = 0;
		iTable[1] = 1;
		iTable[2] = 0;
		iTable[3] = 0;
	}
	else if(cC == 'g' || cC == 'G'){
		iTable[0] = 0;
		iTable[1] = 0;
		iTable[2] = 1;
		iTable[3] = 0;
	}
	else if(cC == 'c' || cC == 'C'){
		iTable[0] = 0;
		iTable[1] = 0;
		iTable[2] = 0;
		iTable[3] = 1;
	}
	else if(cC == 'b' || cC == 'B'){  /* b = g|c|t */
		iTable[0] = 0;
		iTable[1] = 1;
		iTable[2] = 1;
		iTable[3] = 1;
	}
	else if(cC == 'd' || cC == 'D'){  /* d = g|a|t */
		iTable[0] = 1;
		iTable[1] = 1;
		iTable[2] = 1;
		iTable[3] = 0;
	}
	else if(cC == 'h' || cC == 'H'){  /* h = t|c|a */
		iTable[0] =1;
		iTable[1] =1;
		iTable[2] =0;
		iTable[3] =1;
	}
	else if(cC == 'v' || cC == 'V'){  /* v = g|a|c */
		iTable[0] = 1;
		iTable[1] = 0;
		iTable[2] = 1;
		iTable[3] = 1;
	}
	else if(cC == 'y' || cC == 'Y'){  /* y = t|c */
		iTable[0] = 0;
		iTable[1] = 1;
		iTable[2] = 0;
		iTable[3] = 1;
	}
	else if(cC == 'r' || cC == 'R'){  /* r = a|g */
		iTable[0] = 1;
		iTable[1] = 0;
		iTable[2] = 1;
		iTable[3] = 0;
	}
	else if(cC == 'm' || cC == 'M'){  /* m = a|c */
		iTable[0] = 1;
		iTable[1] = 0;
		iTable[2] = 0;
		iTable[3] = 1;
	}
	else if(cC == 'k' || cC == 'K'){  /* k = g|t */
		iTable[0] = 0;
		iTable[1] = 1;
		iTable[2] = 1;
		iTable[3] = 0;
	}
	else if(cC == 'w' || cC == 'W'){  /* w = a|t */
		iTable[0] = 1;
		iTable[1] = 1;
		iTable[2] = 0;
		iTable[3] = 0;
	}
	else if(cC == 's' || cC == 'S'){  /* s = g|c */
		iTable[0] = 0;
		iTable[1] = 0;
		iTable[2] = 1;
		iTable[3] = 1;
	}
	else if(cC == 'n' || cC == 'N' || cC == '-'){
		iTable[0] = 1;
		iTable[1] = 1;
		iTable[2] = 1;
		iTable[3] = 1;
	}
	else{
		cerr << "Error: unrecognized sequence character '" << cC <<"'.\n";
		exit(-1);
	}
	return iTable;
}

string Utilities::checkForConflicts(string & filename){
	ifstream testForExistence;
	string tempFilename = filename;
	unsigned int index = 1;

	testForExistence.open(tempFilename.c_str());
	testForExistence.close();

	while(!testForExistence.fail()){
		if(index > 100){
			//prevents infinite loops
			cerr << "Error: could not find suitable non-existent filename for filename '" << filename <<"'\n";
			exit(-1);
		}
		stringstream ss;
		ss << filename + "." << index;
		index++;
		tempFilename = ss.str();

		testForExistence.clear();
		testForExistence.open(tempFilename.c_str());
		testForExistence.close();
	}

	testForExistence.clear();
	return tempFilename;

}

char Utilities::reverseInterpret(int i){
	assert(i >= 0);
	assert(i < 4);
	string table = "AUGC";
	//cout << i << ":" << table.at(i);
	return table.at(i);
}

bool Utilities::validBP(char a, char b){
	//note: there are probably quicker ways to implement this

	b = toupper(b);
	a = toupper(a);

	//note: this is no longer necessary, since conversion is done when
	//sequences are read in.
	if (a == 'T'){
		a = 'U';
	}
	else if(b == 'T'){
		b = 'U';
	}

	string temp("");
	temp += a;
	temp += b;

	//note: compare returns 0 when the strings are equal, hence the !
	if (!temp.compare("AU") || !temp.compare("UA") ){
		return true;
	}

	if (!temp.compare("GU") || !temp.compare("UG")){
		return true;
	}

	if (!temp.compare("GC") || !temp.compare("CG")){
		return true;
	}

	return false;
}

void Utilities::TrimSpaces(string& str)
{
	//code from: http://sarathc.wordpress.com/2007/01/31/how-to-trim-leading-or-trailing-spaces-of-string-in-c/
	//changes: now trims endlines as well as spaces

	// Trim Both leading and trailing spaces
	size_t startpos = str.find_first_not_of(" \t\n"); // Find the first character position after excluding leading blank spaces
	size_t endpos = str.find_last_not_of(" \t\n"); // Find the first character position from reverse af

	// if all spaces or empty return an empty string
	if(( string::npos == startpos ) || ( string::npos == endpos))
	{
		str = "";
	}
	else
		str = str.substr( startpos, endpos-startpos+1 );
	/*
	     // Code for  Trim Leading Spaces only
	     size_t startpos = str.find_first_not_of(� \t�); // Find the first character position after excluding leading blank spaces
	     if( string::npos != startpos )
	         str = str.substr( startpos );
	 */

	/*
	     // Code for Trim trailing Spaces only
	     size_t endpos = str.find_last_not_of(� \t�); // Find the first character position from reverse af
	     if( string::npos != endpos )
	         str = str.substr( 0, endpos+1 );
	 */

}
