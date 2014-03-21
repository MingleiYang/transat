/*
 * Copyright: Nicholas P Wiebe (2008-2009) and Irmtraud M Meyer (2008-2009)
 * License: licensed under the GNU General Public License version 3 (GPLv3)
 * HelixCore.h
 * Class for storing Helix 'cores' - this concept is not used in the most
 * recent strategies.
 */

#ifndef HELIXCORE_H_
#define HELIXCORE_H_

#include <string>
#include <vector>

using namespace std;

class HelixCore {
public:
	HelixCore(int participatingSeqs, int pos5, int pos3);
	virtual ~HelixCore();

	void growOut(); //add bp outside
	void growIn(); //add bp inside

	string dotBracket(int alignmentLength);

    int getPos5() const
    {
        return pos5;
    }

    int getPos3() const
    {
        return pos3;
    }

    int getLength() const
    {
        return length;
    }

	string printAlignmentAtCore(vector<string*> & alignedSeqs);

	int* countCovarying();


private:
	int participatingSeqs;
	//not currently storing which seqs are participating, just
	//the number of seqs have participating sequences

	int pos5; //5'-most base
	int pos3; //pairing partner of 5'-most base
	int length;

	// example:
	// --------(((---)))----
	// ----pos5^---pos3^----
	//length = 3

};

#endif /* HELIXCORE_H_ */
