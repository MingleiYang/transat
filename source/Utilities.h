/*
 * Copyright: Nicholas P Wiebe (2008-2009) and Irmtraud M Meyer (2008-2009)
 * License: licensed under the GNU General Public License version 3 (GPLv3)
 * Utilities.h
 */

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <vector>
#include <string>

using namespace std;

class Utilities {

public:

	static int iCheckLongHelix(const int* a, const int b);
	static bool iCheckLongHelix(vector<int>& Bps, int pos);
	static int iAt(char cC);

	static void ReadStack(double eStack[4][4][4][4]);

	/**
	 * similar to STL lowerBound, except returns an index instead of an
	 * iterator
	 *
	 */
	static int lowerBound(vector<double> & v, int begin, int end, double value);

	static vector<int> interpret(char a); //similar to iAt, but allows full set of nucleotide symbols

	/**
	 * retrieves the corresponding nucleotide symbol for an integer i
	 * The order is as in interpret (and different from iAt)
	 */
	static char reverseInterpret(int i);

	/**
	 * Checks to see if file with given filename exists
	 * If no, returns the filename unchanged
	 * If yes, returns the filename with an integer prefix added.
	 * Returned filename does not currently exist. (Possibly a race condition?)
	 */
	static string checkForConflicts(string & filename);

	/**
	 * Checks if the base pair formed by nucleotides a and b is a
	 * canonical bp (AU, GC, or GU)
	 */
	static bool validBP(char a, char b);

	static void TrimSpaces(string& str);

//no constructor... just static methods
private:
	Utilities();
	virtual ~Utilities();

};

#endif /* UTILITIES_H_ */
