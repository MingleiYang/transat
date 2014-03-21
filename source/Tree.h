/*
 * Copyright: Nicholas P Wiebe (2008-2009) and Irmtraud M Meyer (2008-2009)
 * License: licensed under the GNU General Public License version 3 (GPLv3)
 * Tree.h
 * Node of phylogenetic tree
 * Possibly this should be split into Tree, storing root information, and TreeNode
 */

#ifndef TREE_H_
#define TREE_H_

#include <string>
#include <list>
#include <map>
#include <vector>

using namespace std;

class Alignment;

class Tree {
public:
	Tree(string treeFileName);
	Tree(list<string*> newickStringTokens);
	Tree(int leaves, double length);

	virtual ~Tree();

	Tree* leftChild;
	Tree* rightChild;
	Tree* parent;
	double branchLength;

	string seqName; //for leaves only - otherwise = "";

	double calcFelsDouble(Alignment & a, int pos5, int pos3);
	double calcFelsSingle(Alignment & a, int pos);

	/*
	 * returns the length of the branch + the total length of all subtrees
	 */
	double totalLength();

	map<Tree*, unsigned int> getLeaf2SeqMap(vector<string*> & SeqNames);
	map<Tree*, unsigned int> getLeaf2SeqMap(vector<string> & SeqNames); //ugly hack for alignment generator...
	vector<string> getSeqNames();

	string newickString();


	double matrix16[16][16];
	double matrix4[4][4];
	double fels[16];

	static bool nonGapPair;

private:

	void createChildren(list<string*> & newickStringTokens);
	void setUpMatrices();

	vector<int> interpretNonPairingGap(vector<int> & nongap);

	static int leafCount;
	static bool felsDoubleUnderflow;
	static bool felsSingleUnderflow;

};

#endif /* TREE_H_ */
