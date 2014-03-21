/*
 * Copyright: Nicholas P Wiebe (2008-2009) and Irmtraud M Meyer (2008-2009)
 * License: licensed under the GNU General Public License version 3 (GPLv3)
 */

#include "HelixFinder.h"
#include "Alignment.h"
#include "Utilities.h"
#include "AlignedHelix.h"
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include "ShuffledAlignment.h"
#include <cassert>
#include <utility>
#include <list>

Realigner HelixFinder::realign = NO_REALIGN;
bool HelixFinder::verbose_out = true;

HelixFinder::HelixFinder(Alignment * a) : alignment(a){
}

HelixFinder::~HelixFinder() {

}

void HelixFinder::findAllHelices(){

	//at this point, store helices in bins according to their outer base-pair
	vector<list<AlignedHelix> > helicesByOuterBP;

	//initialize storage container
	for(unsigned int i = 0; i < alignment->alignedStruct.size() - 1; i++){
		for(unsigned int j = i + 1; j < alignment->alignedStruct.size(); j++){
			helicesByOuterBP.push_back(list<AlignedHelix>());
			assert(convert(i,j) == (int)helicesByOuterBP.size()- 1);
		}
	}

	//find helices in all sequences
	for(unsigned int i = 0; i < alignment->alignedSeqs.size(); i++){
		findAllHelices(i, helicesByOuterBP);
	}

	//store all helices found:
	for(unsigned int j = 0; j < helicesByOuterBP.size(); j++){
		for(list<AlignedHelix>::iterator i = helicesByOuterBP[j].begin(); i != helicesByOuterBP[j].end(); i++){
			//cout << i->dotBracket(alignment->alignedStruct.size()) << "\t" << j << endl;
			helices.push_back(*i);
		}
	}

}

void HelixFinder::findAllHelices(int seqIndex, vector<list<AlignedHelix> > & helicesByOuterBP){

	int **DP;
	int i, j;

	//double eStack[4][4][4][4];
	//Utilities::ReadStack(eStack);

	string Seq = *(alignment->seqs[seqIndex]);
	int iL = Seq.length();

	vector<int> P = *(alignment->seqStructs[seqIndex]);
	// Initialize values

	// allocate and initialize memory for dynamic programming

	DP = (int**) malloc(sizeof(int*) * iL);
	for (i = 0; i < iL; i++) {
		DP[i] = (int*) malloc(sizeof(int) * iL);
	}

	for (i = 0; i < iL; i++) {
		for (j = 0; j < iL; j++) {
			DP[i][j] = 0;
		}
	}

	// dynamic programming

	for (j = MIN_DIST; j < iL; j++) {
		for (i = j - MIN_DIST; i >= 0; i--) {

			assert(i>= 0 && i < iL-MIN_DIST+1);

			if ((Seq[i] == 'A' && Seq[j] == 'U') || (Seq[i] == 'U' && Seq[j] == 'A') || (Seq[i] == 'C' && Seq[j] == 'G') || (Seq[i] == 'G' && Seq[j] == 'C') || (Seq[i] == 'G' && Seq[j] == 'U') || (Seq[i] == 'U' && Seq[j] == 'G')) {
				assert(i+1 < iL-MIN_DIST+1);
				DP[i][j] = DP[i + 1][j - 1] + 1;
			} else {
				DP[i][j] = 0;
			}

			if ((DP[i][j] == 0) && (DP[i + 1][j - 1] > Alignment::minStemLength)) {

				//found a helix:
				Helix h(i+1, j-1, DP[i+1][j-1]);
				int loc = convert(alignment->seq2AlignmentMap[seqIndex]->at(i+1), alignment->seq2AlignmentMap[seqIndex]->at(j-1));

				//try to insert it into existing helices;
				bool found = false;
				for(list<AlignedHelix>::iterator it = helicesByOuterBP[loc].begin(); it != helicesByOuterBP[loc].end(); it++){
					found = it->insert(h, seqIndex, *alignment->seq2AlignmentMap[seqIndex]);
					if(found){
						break;
					}
				}

				//not found -> create new AlignedHelix and add it in
				if(!found){
					AlignedHelix alignedH(h, seqIndex, *alignment->seq2AlignmentMap[seqIndex]);
					helicesByOuterBP[loc].push_back(alignedH);
				}

			}
		}
	}

	//find helices on edges of matrix
	for(j = 0; j < iL; j++){

		if (DP[0][j] > Alignment::minStemLength){
			//found a helix:
			Helix h(0, j, DP[0][j]);
			int loc = convert(alignment->seq2AlignmentMap[seqIndex]->at(0), alignment->seq2AlignmentMap[seqIndex]->at(j));

			//try to insert it into existing helices;
			bool found = false;
			for(list<AlignedHelix>::iterator it = helicesByOuterBP[loc].begin(); it != helicesByOuterBP[loc].end(); it++){
				found = it->insert(h, seqIndex, *alignment->seq2AlignmentMap[seqIndex]);
				if(found){
					break;
				}
			}

			//not found -> create new AlignedHelix and add it in
			if(!found){
				AlignedHelix alignedH(h, seqIndex, *alignment->seq2AlignmentMap[seqIndex]);
				helicesByOuterBP[loc].push_back(alignedH);
			}

		}

		if (DP[j][iL-1] > Alignment::minStemLength && j != 0){ //the j==0 case is caught by the above condition (when j = iL-1)
			//found a helix:
			Helix h(j, iL-1, DP[j][iL-1]);
			int loc = convert(alignment->seq2AlignmentMap[seqIndex]->at(j), alignment->seq2AlignmentMap[seqIndex]->at(iL-1));

			//try to insert it into existing helices;
			bool found = false;
			for(list<AlignedHelix>::iterator it = helicesByOuterBP[loc].begin(); it != helicesByOuterBP[loc].end(); it++){
				found = it->insert(h, seqIndex, *alignment->seq2AlignmentMap[seqIndex]);
				if(found){
					break;
				}
			}

			//not found -> create new AlignedHelix and add it in
			if(!found){
				AlignedHelix alignedH(h, seqIndex, *alignment->seq2AlignmentMap[seqIndex]);
				helicesByOuterBP[loc].push_back(alignedH);
			}
		}

	}

	for (i = 0; i < iL; i++) {
		free(DP[i]);
	}
	free(DP);

	return;
}

int HelixFinder::convert(int pos5, int pos3){
	return pos5 * alignment->alignedStruct.size() + pos3 - ((pos5+2)*(pos5+1))/2;
}

void HelixFinder::allHelicesPvalueTable(int randomSamples, Tree & tree, bool doPvalues){

	assert(randomSamples > 0);
	cout << "Pvalue\tLogLikeRatio\tPairedLogLikelihood\tUnpairedLogLikelihood\tMeanFracCanonicalBP\tCovariance\tConservation\tLength\tMidpoint\tNormalizedMidpoint\tAppearances\tIsConsensusHelix\tIsPartialConsensusHelix\tIsCompetingHelix\tConsensusBps\tStructure"
	<< "\tAlignmentSize\tAlignmentLength\tTreeLength\tBps"
	<< "\tCis\tTrans"
	<< "\tCis5\tCis3\tTrans5\tTrans3\tMid5\tMid3"
	<< "\tnewPvalue"
	<< endl;

	if(helices.empty()){
		cerr << "No helices found - Pvalue Table empty\n";
		return;
	}

	//initialize pvalues
	vector<double> pvalues(helices.size(), 0.0);
	vector<double> newPvalues(helices.size(), 0.0);
	vector<unsigned int> lowerHelixCount(helices.size(), 0);
	unsigned int nullHelixCount = 0;
	vector<double> randLogLikes;
	//if we don't do p value calculation, p-value column will be all zeros

	Alignment * startingAlignment;

	if(realign == TCOFFEE){
		vector<string> realigned_seqs = alignment->realignTcoffee(tree);
		vector<int> emptyStruct(realigned_seqs.begin()->length(),-1);
		vector<string> names;
		//TODO: it would really be nice to get rid of the vectors of pointers in Alignment class...
		for(unsigned int i = 0; i < alignment->seqNames.size(); i++){
			names.push_back(*alignment->seqNames[i]);
		}
		startingAlignment = new Alignment(names, realigned_seqs, emptyStruct);
	}
	else{
		startingAlignment = alignment;
	}

	if(doPvalues){
		for(int i = 0; i < randomSamples; i++){


			ShuffledAlignment randomizedAlignment(*startingAlignment);
			HelixFinder randomHelices(&randomizedAlignment);

			randomHelices.findAllHelices();
			for(unsigned int j = 0; j < randomHelices.helices.size(); j++){
				randLogLikes.push_back(randomHelices.helices[j].logLikeRatio(randomizedAlignment, tree));
			}

			sort(randLogLikes.begin(), randLogLikes.end());

			if(!randLogLikes.empty()){
				int lowerbound;
				for(unsigned int j = 0; j < pvalues.size(); j++){
					lowerbound = Utilities::lowerBound(randLogLikes, 0, randLogLikes.size(), helices[j].logLikeRatio(*alignment, tree));
					lowerHelixCount[j] += lowerbound;
					pvalues[j] += 1 - lowerbound/(double) randLogLikes.size();
				}
			}
			//if empty, implicitly add 0 to each pvalue

			nullHelixCount += randLogLikes.size();

			randLogLikes.clear();
		}

		//take average over all samples
		assert(nullHelixCount >0);
		for(unsigned int j = 0 ; j < pvalues.size(); j++){
			newPvalues[j] = 1 - lowerHelixCount[j] / (double)nullHelixCount;
			pvalues[j] = pvalues[j] / randomSamples;
		}
	}

	double treeLength = tree.totalLength();

	for(unsigned int j = 0; j < helices.size(); j++){
//		cout << helices[j].isConsensusHelix(*alignment) << "\t";
		cout << pvalues[j] << "\t";
		cout << helices[j].logLikeRatio(*alignment, tree) << "\t";
		cout << helices[j].pairedLikelihood(*alignment, tree) << "\t";
		cout << helices[j].unpairdLikelihood(*alignment, tree) << "\t";
		cout << helices[j].canonicalBP(*alignment) << "\t";
		cout << helices[j].covariance(*alignment) << "\t";
		cout << helices[j].conservation(*alignment) << "\t";
		cout << helices[j].length() << "\t";
		cout << helices[j].midpoint() << "\t";
		cout << helices[j].midpoint(*alignment) << "\t";
		cout << helices[j].appearances() << "\t";
		cout << helices[j].isConsensusHelix(*alignment) << "\t";
		cout << helices[j].isPartialConsensusHelix(*alignment) << "\t";
		cout << helices[j].isCompetingHelix(*alignment) << "\t";
		cout << helices[j].consensusBps(*alignment) << "\t";
		if(verbose_out){
			cout << helices[j].dotBracket(alignment->alignedStruct.size());
		}
		else{
			cout << 0;
		}
		cout << "\t" << alignment->seqs.size();
		cout << "\t" << alignment->alignedStruct.size();
		cout << "\t" << treeLength;
		cout << "\t" << helices[j].bpsString();

		pair<double, double> cisTrans = helices[j].cisTransScore(*alignment);
		cout << "\t" << cisTrans.first;
		cout << "\t" << cisTrans.second;

		double cis5, cis3, trans5, trans3, mid5, mid3;
		helices[j].competeScore(*alignment, cis5, cis3, trans5, trans3, mid5, mid3);
		cout << "\t" << cis5;
		cout << "\t" << cis3;
		cout << "\t" << trans5;
		cout << "\t" << trans3;
		cout << "\t" << mid5;
		cout << "\t" << mid3;

		cout << "\t" << newPvalues[j];

		cout << endl;
	}

}

void HelixFinder::balancedSparseHelixTable(Tree & tree){
	vector<AlignedHelix> fakeHelices;
	for(unsigned int i = 0; i < helices.size(); i++){
		if(helices[i].isConsensusHelix(*alignment)){
			cout <<  helices[i].isConsensusHelix(*alignment) << " ";
			cout << "1:" << helices[i].pairedLikelihood(*alignment, tree) << " ";
			cout << "2:" << helices[i].unpairdLikelihood(*alignment, tree) << " ";
			cout << "3:" << helices[i].canonicalBP(*alignment) << " ";
			cout << "4:" << helices[i].covariance(*alignment) << " ";
			cout << "5:" << helices[i].conservation(*alignment);
			cout << endl;
		}
		else{
			fakeHelices.push_back(helices[i]);
		}
	}
	srand((unsigned)time(0));
	for(unsigned int i = 0; i < helices.size() - fakeHelices.size(); i++){
		AlignedHelix temp = fakeHelices[i];
		unsigned int random = rand();

		random = i + random % (fakeHelices.size()-i);

		fakeHelices[i] = fakeHelices[random];
		fakeHelices[random] = temp;
		cout <<  fakeHelices[i].isConsensusHelix(*alignment) << " ";
		cout << "1:" << fakeHelices[i].pairedLikelihood(*alignment, tree) << " ";
		cout << "2:" << fakeHelices[i].unpairdLikelihood(*alignment, tree) << " ";
		cout << "3:" << fakeHelices[i].canonicalBP(*alignment) << " ";
		cout << "4:" << fakeHelices[i].covariance(*alignment) << " ";
		cout << "5:" << fakeHelices[i].conservation(*alignment);
		cout << endl;
	}
}

void HelixFinder::sparseHelixTable(Tree & tree){
	for(unsigned int i = 0; i < helices.size(); i++){

		cout <<  helices[i].isConsensusHelix(*alignment) << " ";
		cout << "1:" << helices[i].pairedLikelihood(*alignment, tree) << " ";
		cout << "2:" << helices[i].unpairdLikelihood(*alignment, tree) << " ";
		cout << "3:" << helices[i].canonicalBP(*alignment) << " ";
		cout << "4:" << helices[i].covariance(*alignment) << " ";
		cout << "5:" << helices[i].conservation(*alignment);
		cout << endl;

	}
}

void HelixFinder::findTrueHelices(){

}

void HelixFinder::coverage(double & coverage, double & exactMatchCoverage, double trueHelixCutoff){
	assert(trueHelixCutoff >=0);
	assert(trueHelixCutoff <=1);

	vector<bool> covered(alignment->alignedStruct.size(), false);
	vector<bool> exactCovered(alignment->alignedStruct.size(), false);


	for(unsigned int i = 0 ; i < helices.size(); i++){
		if(helices[i].isConsensusHelix(*alignment)){
			for(unsigned int j = 0; j < helices[i].bps.size(); j++){
				exactCovered[helices[i].bps[j].first] = true;
				exactCovered[helices[i].bps[j].second] = true;
			}
		}
		if(helices[i].consensusBps(*alignment) / (double) helices[i].length() > trueHelixCutoff || (trueHelixCutoff == 1.0 && helices[i].consensusBps(*alignment) == (int)helices[i].length()) ){
			for(unsigned int j = 0; j < helices[i].bps.size(); j++){
				if(alignment->alignedStruct[helices[i].bps[j].first] == helices[i].bps[j].second){
					covered[helices[i].bps[j].first] = true;
					covered[helices[i].bps[j].second] = true;
				}
			}
		}

	}

	int exactMatches = 0;
	int matches = 0;
	int pairedBaseCount = 0;

	for(unsigned int i = 0; i < covered.size(); i++){
		if(alignment->alignedStruct[i] >= 0){
			pairedBaseCount++;
			if(exactCovered[i]){
				assert(covered[i]);
				exactMatches++;
				matches++;
			}
			else if(covered[i]){
				matches++;
			}
		}
	}
	coverage = matches / (double)pairedBaseCount;
	exactMatchCoverage = exactMatches / (double) pairedBaseCount;

}

void HelixFinder::findAllHelicesGrow(Tree & tree, double seedThreshold, double growThreshold){

	//matrix that stores which bps score high enough to be considered part of a helix
	//partOfHelix[i][j-i] gives the value for bp i-j (i is 5' of j)
	int alignmentLength = alignment->alignedStruct.size();

	vector<vector<bool> > partOfHelix(alignmentLength);
	for(unsigned int i = 0; i < partOfHelix.size(); i++){
		partOfHelix[i] = vector<bool>(partOfHelix.size() - i, false);
	}

	//normalizing factor
	double treeLength = tree.totalLength();

//	list<pair<int, int> > outerGrow;
//	list<pair<int, int> > innerGrow;
	list<pair<int, int> > growPoints;

	//find seed points
	for(int i = 0; i < alignmentLength-MIN_DIST; i++){
		for(int j = i + MIN_DIST; j < alignmentLength; j++){
			if(alignment->logLikelihood(i,j, tree) / treeLength > seedThreshold){
				partOfHelix[i][j-i] = true;
//				if(i > 0 && j+1 < alignment->alignedStruct.size()){
//					outerGrow.push_back(pair<int,int>(i,j));
//				}
//				if(i+MIN_DIST < j){
//					innerGrow.push_back(pair<int,int>(i,j));
//				}
				growPoints.push_back(pair<int,int>(i,j));
			}
		}
	}

	//grow from seed points
	pair<int, int> candidate;
	pair<int, int> current;

	while(!growPoints.empty()){
		current = growPoints.front();
		//outer bp growth:
		if(current.first > 0 && current.second + 1 < alignmentLength){
			candidate = pair<int,int>(current.first -1, current.second +1);
			if(!partOfHelix[candidate.first][candidate.second - candidate.first] && alignment->logLikelihood(candidate.first, candidate.second , tree) / treeLength > growThreshold){
				partOfHelix[candidate.first][candidate.second - candidate.first] = true;
				growPoints.push_back(candidate);
			}
		}


		if(current.first + 1 + MIN_DIST <= current.second - 1){
			candidate = pair<int,int>(current.first +1, current.second -1);
			if(!partOfHelix[candidate.first][candidate.second - candidate.first] && alignment->logLikelihood(candidate.first, candidate.second , tree) / treeLength > growThreshold){
				partOfHelix[candidate.first][candidate.second - candidate.first] = true;
				growPoints.push_back(candidate);
			}
		}
		growPoints.pop_front();
	}

	for(int i = 0; i < alignmentLength; i++){
		for(int k = 0; k < i; k++){
			cout << "\t";
		}
		for(int j = i; j < alignmentLength-1; j++){
			cout << partOfHelix[i][j-i] << "\t";
		}
		cout << partOfHelix[i][alignmentLength-1-i];
		cout << "\n";
	}

	int helixLength = 0;
	for(int i = 0; i < alignmentLength; i++){
		int j;
		for(j = 0; i-j >= 0 && i+j < alignmentLength; j++){
			if(partOfHelix[i-j][2*j]){
				helixLength++;
			}
			else{
				if(helixLength > Alignment::minStemLength){
					AlignedHelix helix(i-j+1, i+j-1,helixLength);
					helices.push_back(helix);
				}
				helixLength = 0;
			}
		}
		if(helixLength > Alignment::minStemLength){
			AlignedHelix helix(i-j+1, i+j-1,helixLength);
			helices.push_back(helix);
		}

		helixLength = 0;
	}

}
