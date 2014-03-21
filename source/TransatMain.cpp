/*
 * Copyright: Nicholas P Wiebe (2008-2009) and Irmtraud M Meyer (2008-2009)
 * License: licensed under the GNU General Public License version 3 (GPLv3)
 * Main function, runs Transat analysis on all aligned helices.
 */

#include <iostream>
#include "Alignment.h"
#include <string>
#include "Utilities.h"
#include <fstream>
#include "Tree.h"
#include <cassert>
#include <math.h>
#include "AlignmentGenerator.h"
#include "HelixFinder.h"
#include "AlignedHelix.h"
#include "ShuffledAlignment.h"
#include <algorithm>

using namespace std;

int main(int argc, char** argv) {



	string treeFile = "";
	bool coverageInfo = false;
	bool useTree = false;
	int randomTrials = 500;
	bool pVals = true;
	bool bpTable = false;
	bool growHelices = false;
	string filename = "";
	string structFilename = "";
	bool noStruct = false;

	for(int i = 1; i < argc; i++){
		string temp = argv[i];
		//cout << "arg: " << temp << endl;
		if(temp.compare("-fasta") == 0){
			i++;
			if(i < argc){
				filename = argv[i];
			}
			else{
				cerr << "Error: expecting filename following \"-fasta\" argument\n";
				exit(-1);
			}
		}
		else if (temp.compare("-ct") == 0){
			i++;
			if(i < argc){
				structFilename = argv[i];
			}
			else{
				cerr << "Error: expecting filename following \"-ct\" argument\n";
				exit(-1);
			}
		}
		else if (temp.compare("-minSL") == 0){
			i++;
			if(i < argc){
				Alignment::minStemLength = atoi(argv[i]);
			}
			else{
				cerr << "Error: expecting integer following \"-minSL\" argument\n";
				exit(-1);
			}
		}
		else if (temp.compare("-tree") == 0){
			i++;
			if(i < argc){
				treeFile = argv[i];
				useTree = true;
			}
			else{
				cerr << "Error: expecting filename following \"-tree\" argument\n";
				exit(-1);
			}
		}
		else if (temp.compare("-randomize") == 0){
			i++;
			if(i < argc){
				randomTrials = atoi(argv[i]);
			}
			else{
				cerr << "Error: expecting integer following \"-randomize\" argument\n";
				exit(-1);
			}
		}
		else if(temp.compare("-noPvalues") == 0){
			//if this argument is given, no null distributions will be generated, and so no
			//pvalues will be calculated. The Pvalue column in the output will be filled with zeros
			pVals = false;
		}
		else if(temp.compare("-bpTable") == 0){
			bpTable = true;
		}
		else if(temp.compare("-realign") == 0){
			HelixFinder::realign = TCOFFEE;
		}
		else if(temp.compare("-coverage") == 0){
			coverageInfo = true;
		}
		else if(temp.compare("-grow") == 0){
			growHelices = true;
		}
		else if(temp.compare("-nonGapPair") == 0){
			Tree::nonGapPair = !Tree::nonGapPair;
		}
		else if(temp.compare("-noDB") == 0){
			//no dot-bracket
			HelixFinder::verbose_out = false;
		}
		else if(temp.compare("-noStruct") == 0){
			//no structure included
			noStruct = true;
		}
		else{
			cerr << "Error: unrecognized argument '" << temp <<"'\n";
			exit(-1);
		}
	}

	//check for missing args
	bool missingArgs = false;
	if(filename.compare("") == 0){
		cerr << "Error: must supply an alignment file (use '-fasta <filename>')\n.";
		missingArgs = true;
	}
	if(treeFile.compare("") == 0){
		cerr << "Error: must supply a tree file (use '-tree <filename>')\n";
		missingArgs = true;
	}


	if(missingArgs){
		exit(-1);
	}


	Alignment * a;
	if (!noStruct && structFilename.compare("") == 0){
		a = new Alignment(filename);
	}
	else{
		a = new Alignment(filename, structFilename);
	}



	assert(useTree);
	Tree root(treeFile);
	//sanity check: make sure all sequences are in tree:
	map<Tree*, unsigned int> leaf2SeqMap = root.getLeaf2SeqMap(a->seqNames);
	assert(leaf2SeqMap.size() == a->seqNames.size());
	if(bpTable){
		a->sparseBpTable(root);
	}
	else if(coverageInfo){
		HelixFinder hf(a);
		hf.findAllHelices();

		double coverage;
		double exactCoverage;
		hf.coverage(coverage, exactCoverage, 0.5);
		cout << coverage << "\t" << exactCoverage << endl;

	}
	else{
		HelixFinder hf(a);
		if(growHelices){
			hf.findAllHelicesGrow(root, 0, -5);
		}
		else{
			hf.findAllHelices();
		}
		hf.allHelicesPvalueTable(randomTrials, root, pVals);
	}

	cerr << "done!\n";
	delete a;
	return 0;
}
