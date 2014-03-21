/*
 * Copyright: Nicholas P Wiebe (2008-2009) and Irmtraud M Meyer (2008-2009)
 * License: licensed under the GNU General Public License version 3 (GPLv3)
 * Alignment.cpp
 */

#include "Alignment.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cassert>
#include "Utilities.h"
#include <math.h>
#include <algorithm>
#include "Tree.h"
#include <stack>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include "ShuffledAlignment.h"

string Alignment::TEMP_FASTA_FILENAME = "temp.fasta";
string Alignment::TEMP_OUTPUT_FASTA_FILENAME = "temp.out.fasta";
const string Alignment::REALIGNER_LOCATION = "Realigner/bin";
const string Alignment::TCOFFEE_LOC = "t_coffee";

//string Alignment::JAVA_LOC = "/cs/local/lib/pkg/jdk-1.6.0.02/bin/java";
string Alignment::JAVA_LOC = "java";


int Alignment::minStemLength = 8;
bool Alignment::printHeaders = true;

Alignment::Alignment() {
	consensusBP = NULL;
	gapFraction = NULL;
	seqCons = NULL;

}

Alignment::Alignment(vector<string> & names, vector<string> & alignment, vector<int> & structure){
	//copy names
	for(unsigned int i = 0; i < names.size(); i++){
		seqNames.push_back(new string(names[i]));
	}

	//copy alignment:
	for(unsigned int i = 0; i < names.size(); i++){
		alignedSeqs.push_back(new string(alignment[i]));
	}

	//copy structure:
	alignedStruct = structure;

	//extract unaligned seqs and make map of unaligned seq positions to aligned seq positions
	vector<int>* temp;
	vector<int>* temp2;
	string * seq;

	unsigned int alignedSeqLength = alignedSeqs[0]->length();

	for (unsigned int i = 0; i < alignedSeqs.size(); i++){

		assert(alignedSeqLength == alignedSeqs[i]->length()); //sanity check
		temp = new vector<int>();
		temp2 = new vector<int>();
		seq = new string();

		for (unsigned int j = 0; j < alignedSeqLength; j++){
			if (alignedSeqs[i]->at(j) != '-'){
				temp->push_back(j);
				*seq = *seq + alignedSeqs[i]->at(j);
				temp2->push_back(temp->size() - 1);
			}
			else{
				temp2->push_back(-1);
			}
		}
		seqs.push_back(seq);
		seq2AlignmentMap.push_back(temp);
		alignment2SeqMap.push_back(temp2);

	}

	fillSeqStructs();

	labelHelices2();
	initializeStatsMatrix();

	consensusBP = NULL;
	gapFraction = NULL;
	seqCons = NULL;
}

Alignment::Alignment(string& filename) {
	readAlignment(filename, true);
	labelHelices2();
	initializeStatsMatrix();

	consensusBP = NULL;
	gapFraction = NULL;
	seqCons = NULL;

}

Alignment::Alignment(string & filename, string & structFilename)
{
	readAlignment(filename);
	if(structFilename.compare("")){
		readStruct(structFilename);
	}
	else{
		//no struct file
		emptyStruct();
	}
	labelHelices2();
	initializeStatsMatrix();

	consensusBP = NULL;
	gapFraction = NULL;
	seqCons = NULL;
}

void Alignment::initializeStatsMatrix()
{
	//delete StatsMatrix
	for (unsigned int i = 0; i < StatsMatrix.size(); i++){
		for (unsigned int j = 0; j < StatsMatrix[i]->size(); j++){
			for(unsigned int k = 0; k < StatsMatrix[i]->at(j)->size(); k++){
				delete StatsMatrix[i]->at(j)->at(k);
			}
			StatsMatrix[i]->at(j)->clear();
			delete StatsMatrix[i]->at(j);
		}
		StatsMatrix[i]->clear();
		delete StatsMatrix[i];
	}
	StatsMatrix.clear();

	//delete vector<vector<CompetingHelix*>* > competingHelices
	for (unsigned int i = 0; i < competingHelices.size(); i++){
		for (unsigned int j = 0; j < competingHelices[i]->size(); j++){
			delete competingHelices[i]->at(j);
		}
		competingHelices[i]->clear();
		delete competingHelices[i];
	}
	competingHelices.clear();

	//delete vector<vector<double>* > FelsDoubles
	for(unsigned int i = 0; i < felsDoubles.size(); i++){
		delete felsDoubles[i];
	}
	felsDoubles.clear();
	felsSingles.clear();

	//initialize stats matrix
	vector<vector<StatsWrapper*>* >* temp;
	for (unsigned int i = 0; i < seqs.size(); i++){
		temp = new vector<vector<StatsWrapper*>* >();
		for(int j = 0; j < helixNumber; j++){
			temp->push_back(new vector<StatsWrapper*>());
		}
		StatsMatrix.push_back(temp);

		//create one (empty) vector of competing helices for each sequence
		competingHelices.push_back(new vector<CompetingHelix*>());
	}

	//initialize log likelihood tables:
	for(unsigned int i = 0; i < alignedStruct.size(); i++){
		felsSingles.push_back(1);
		vector<double> * temp_vec = new vector<double>();
		//should be an upper-triangular matrix, but isn't currently
		for(unsigned int j = 0; j < alignedStruct.size(); j++){
			temp_vec->push_back(1);
		}
		felsDoubles.push_back(temp_vec);
	}
}


void Alignment::clearAll(){

	//memory management stuff... boring!
	unsigned int i;
	for (i = 0; i < seqNames.size(); i++){
		delete seqNames[i];
	}
	seqNames.clear();

	for (i = 0; i < seqs.size(); i++){
		delete seqs[i];
	}
	seqs.clear();

	for (i = 0; i < alignedSeqs.size(); i++){
		delete alignedSeqs[i];
	}
	alignedSeqs.clear();

	for (i = 0; i < seq2AlignmentMap.size(); i++){
		delete seq2AlignmentMap[i];
	}
	seq2AlignmentMap.clear();

	for (i = 0; i < alignment2SeqMap.size(); i++){
		delete alignment2SeqMap[i];
	}
	alignment2SeqMap.clear();

	for(i = 0; i < seqStructs.size(); i++){
		delete seqStructs[i];
	}
	seqStructs.clear();

	//delete helixLabels
	for(i = 0; i < helixLabels.size(); i++){
		delete helixLabels[i];
	}
	helixLabels.clear();

	//delete StatsMatrix
	for (i = 0; i < StatsMatrix.size(); i++){
		for (unsigned int j = 0; j < StatsMatrix[i]->size(); j++){
			for(unsigned int k = 0; k < StatsMatrix[i]->at(j)->size(); k++){
				delete StatsMatrix[i]->at(j)->at(k);
			}
			StatsMatrix[i]->at(j)->clear();
			delete StatsMatrix[i]->at(j);
		}
		StatsMatrix[i]->clear();
		delete StatsMatrix[i];
	}
	StatsMatrix.clear();

	//delete vector<vector<CompetingHelix*>* > competingHelices
	for (i = 0; i < competingHelices.size(); i++){
		for (unsigned int j = 0; j < competingHelices[i]->size(); j++){
			delete competingHelices[i]->at(j);
		}
		competingHelices[i]->clear();
		delete competingHelices[i];
	}
	competingHelices.clear();

	for(unsigned int i = 0; i < interestingRegions.size(); i++){
		delete interestingRegions[i];
	}
	interestingRegions.clear();

	for(unsigned int i = 0; i < cores.size(); i++){
		for(unsigned int j = 0; j < cores[i]->size(); j++){
			delete cores[i]->at(j);
		}
		delete cores[i];
	}
	cores.clear();

	for(unsigned int i = 0; i < trueHelices.size(); i++){
		delete trueHelices[i];
	}
	trueHelices.clear();

	//delete vector<vector<double>* > FelsDoubles
	for(unsigned int i = 0; i < felsDoubles.size(); i++){
		delete felsDoubles[i];
	}
	felsDoubles.clear();

	delete consensusBP;
	consensusBP = NULL;

	delete gapFraction;
	gapFraction = NULL;

	delete seqCons;
	seqCons = NULL;

}

Alignment::~Alignment() {
	clearAll();
}

string Alignment::seqString()
{
	string out = "";
	for(unsigned int i = 0; i < seqs.size(); i++){
		out += *(seqs[i]);
		out += "\n";
	}
	return out;
}

string Alignment::alignmentString()
{
	//cout << "Alignment Size: " << alignedSeqs.size() << endl;
	string out = "";
	for(unsigned int i = 0; i < alignedSeqs.size(); i++){
		//cout << *(alignedSeqs[i]) << endl;
		out += *(alignedSeqs[i]);
		out += "\n";
	}
	return out;
}


string Alignment::structString()
{
	string out = "";

	int pos;
	for(unsigned int i = 0; i < seqStructs.size(); i++){
		pos = 0;
		for(unsigned int j = 0; j < seqStructs[i]->size(); j++){
			assert(seqStructs[i]->at(j) != pos);
			if(seqStructs[i]->at(j) == -1){

				out += ".";
			}
			else if (seqStructs[i]->at(j) > pos){
				out += "(";
			}
			else{
				out += ")";
			}

			pos++;
		}
		out += "\n";
	}
	return out;
}

string Alignment::gappedStructString()
{
	string out = "";

	int pos;

	assert(seqStructs.size() == alignedSeqs.size());

	for(unsigned int i = 0; i < seqStructs.size(); i++){
		//alignedStruct
		pos = 0;

		for(unsigned int j = 0; j < alignedSeqs[i]->size(); j++){
			if(alignedSeqs[i]->at(j) == '-' ) {
				out += "-";
			} else {
				assert(seqStructs[i]->at(pos) != pos);
				if(seqStructs[i]->at(pos) == -1){
					out += ".";
				}
				else if (seqStructs[i]->at(pos) > pos){
					out += "(";
				}
				else{
					out += ")";
				}
				pos++;
			}


		}
		out += "\n";
	}
	return out;
}


string Alignment::alignedStructString()
{
	string out = "";
	int j = 0;
	for(unsigned int i = 0; i < alignedStruct.size(); i++){
		assert(alignedStruct[i] != j);
		if(alignedStruct[i] == -1){

			out += ".";
		}
		else if (alignedStruct[i] > j){
			out += "(";
		}
		else{
			out += ")";
		}

		j++;
	}
	return out;
}

string Alignment::alignedlabelString(){

	stringstream ss;
	for (unsigned int i = 0; i < alignmentLabels.size(); i++){
		if(alignmentLabels[i] == -1){
			ss << ".";
		}
		else{
			ss << alignmentLabels[i];
		}
	}
	return ss.str();
}

string Alignment::labelString(){


	stringstream ss;

	for (unsigned int i = 0; i < helixLabels.size(); i++){
		for(unsigned int j = 0; j < helixLabels[i]->size(); j++){
			if(helixLabels[i]->at(j) == -1){
				ss << ".";
			}
			else{
				ss << helixLabels[i]->at(j);
			}
		}
		ss << "\n";
	}

	return ss.str();
}

string Alignment::competingHelicesString(){
	stringstream ss;
	for(unsigned int i = 0; i < competingHelices.size(); i++){
		ss << "Sequence " << *seqNames[i] << ": " << competingHelices[i]->size() << " competing helices\n";
		//for (unsigned int j = 0; j < competingHelices[i]->size(); j++){
		//	ss << competingHelices[i]->at(j)->printHelix2(seqs[i]->length()) << "\n";
		//}
	}
	return ss.str();
}

double Alignment::seqConservationAtPosition(int i)
{
	double baseCounts[4];
	char base;
	unsigned int seqIndex;
	double max;

	//set counts to zero
	for(int j = 0; j < 4; j++){
		baseCounts[j] = 0.0;
	}

	for(seqIndex = 0; seqIndex < alignedSeqs.size(); seqIndex++){
		base = alignedSeqs[seqIndex]->at(i);

			//if position is gap, add one to every letter
		//since we are more interested in low conservation regions
		//this is more a more conservative approach
		vector<int> increment = Utilities::interpret(base);
		for (int j = 0; j< 4; j++){
			baseCounts[j] += increment[j];
		}


	}
	//divide by the number of sequences to get
	//probability for each base
	for(int j = 0; j < 4; j++){
		baseCounts[j] = baseCounts[j] / alignedSeqs.size();
	}

	max = 0;
	for (int j = 0; j < 4; j++){
		if (baseCounts[j] > max){
			max = baseCounts[j];
		}
	}

	return max;
}

void Alignment::rankInterestingRegions(){
	for(int i = 0; i < helixNumber; i++){
		findInterestingRegions(i);
	}
	sort(interestingRegions.begin(), interestingRegions.end(), InterestingRegion::comp);

	//for(int i = 0; i < 10 ; i++){
	//	interestingRegions[interestingRegions.size() - i - 1]->print();
	//}

	for(unsigned int i = 0; i < interestingRegions.size(); i++){
		interestingRegions[i]->print();
		//cout << interestingRegions[i]->seqConservation << "\t" << interestingRegions[i]->structConservation << "\n";
	}
}

void Alignment::findInterestingRegions(int trueHelix)
{
	int window_size = 6;

	vector<double> c3perPos;
	vector<double> c5perPos;
	vector<double> t3perPos;
	vector<double> t5perPos;

	vector<double> seqCons;

	//initialize perPos vectors
	for(unsigned int i = 0; i < alignedStruct.size(); i++){
		c3perPos.push_back(0);
		c5perPos.push_back(0);
		t3perPos.push_back(0);
		t5perPos.push_back(0);
		seqCons.push_back(seqConservationAtPosition(i)); //calculate sequence conservation
	}

	//fill perPosVectors... like in competingHelicesHistogramPerHelix
	unsigned int totalSum = 0;

	vector<int> pos;
	vector<int> helixNumCompetingAtPos;
	vector<int> helixNumCompetingAtPos_c3;
	vector<int> helixNumCompetingAtPos_t3;
	vector<int> helixNumCompetingAtPos_c5;
	vector<int> helixNumCompetingAtPos_t5;
	for(unsigned int j = 0; j < numSeqs(); j++) {
		pos.push_back(0);
		helixNumCompetingAtPos.push_back(0);
		helixNumCompetingAtPos_c3.push_back(0);
		helixNumCompetingAtPos_t3.push_back(0);
		helixNumCompetingAtPos_c5.push_back(0);
		helixNumCompetingAtPos_t5.push_back(0);
	}
	for (unsigned int j = 0; j < alignedSeqs[0]->size(); j++){
		for (unsigned int seqIndex = 0; seqIndex < alignedSeqs.size(); seqIndex++){
			helixNumCompetingAtPos[seqIndex] = 0;
			helixNumCompetingAtPos_c3[seqIndex] = 0;
			helixNumCompetingAtPos_t3[seqIndex] = 0;
			helixNumCompetingAtPos_c5[seqIndex] = 0;
			helixNumCompetingAtPos_t5[seqIndex] = 0;
			if ( alignedSeqs[seqIndex]->at(j) != '-' ) {
				int helixLabelatPos = helixLabels[seqIndex]->at(pos[seqIndex]);
				if( trueHelix != helixLabelatPos ) {
					for (int competingSeqIndex = 0; competingSeqIndex < (int)competingHelices[seqIndex]->size(); competingSeqIndex++) {
						//if ( competingHelices[seqIndex]->at(competingSeqIndex)->trueHelixIndex == trueHelix ) {
						if (!StatsMatrix[seqIndex]->at(trueHelix)->at(competingSeqIndex)->isZero() ) {
							assert(helixNumCompetingAtPos[seqIndex] >= 0);
							if(competingHelices[seqIndex]->at(competingSeqIndex)->isPosCompeting(pos[seqIndex])) {
								//helixNumCompetingAtPos[seqIndex] += competingHelices[seqIndex]->at(competingSeqIndex)->isPosCompeting(pos[seqIndex]);
								helixNumCompetingAtPos[seqIndex]++;
								if(StatsMatrix[seqIndex]->at(trueHelix)->at(competingSeqIndex)->cis5() > 0) {
									//cout << "DEBUG cis5" << StatsMatrix[seqIndex]->at(trueHelix)->at(competingSeqIndex)->cis5() << endl;
									helixNumCompetingAtPos_c5[seqIndex]++;
								}
								if(StatsMatrix[seqIndex]->at(trueHelix)->at(competingSeqIndex)->cis3() > 0) {
									helixNumCompetingAtPos_c3[seqIndex]++;
								}
								if(StatsMatrix[seqIndex]->at(trueHelix)->at(competingSeqIndex)->trans5() > 0) {
									helixNumCompetingAtPos_t5[seqIndex]++;
								}
								if(StatsMatrix[seqIndex]->at(trueHelix)->at(competingSeqIndex)->trans3() > 0) {
									helixNumCompetingAtPos_t3[seqIndex]++;
								}
							}
						}
					}
				} else {
					assert(helixNumCompetingAtPos[seqIndex] == 0);
					helixNumCompetingAtPos[seqIndex] = -1;
				}
				pos[seqIndex]++;
			}
		}
		int tmpTrue = 0;
		int tmpComp = 0;
		int c3 = 0;
		int c5 = 0;
		int t3 = 0;
		int t5 = 0;

		for(unsigned int i = 0; i < numSeqs(); i++) {
			if(helixNumCompetingAtPos[i] > 0)
				tmpComp++;
			else if(helixNumCompetingAtPos[i] < 0)
				tmpTrue--;
			if(helixNumCompetingAtPos_c3[i] > 0) { c3++; }
			if(helixNumCompetingAtPos_c5[i] > 0) { c5++; }
			if(helixNumCompetingAtPos_t3[i] > 0) { t3++; }
			if(helixNumCompetingAtPos_t5[i] > 0) { t5++; }
		}
		totalSum += tmpComp;

		c3perPos[j] = (double) c3/numSeqs();
		c5perPos[j] = (double) c5/numSeqs();
		t3perPos[j] = (double) t3/numSeqs();
		t5perPos[j] = (double) t5/numSeqs();

	}


	//check all sequence windows, filling interesting regions vector
	for(int i = 0; i < ((int) alignedStruct.size()) - window_size + 1; i++){
		double avg_cons = 0;
		double avg_c3 = 0;
		double avg_c5 = 0;
		double avg_t3 = 0;
		double avg_t5 = 0;

		for (int j = 0; j < window_size; j++){
			avg_cons += seqCons[i+j];

			avg_c3 += c3perPos[i+j];
			avg_c5 += c5perPos[i+j];
			avg_t3 += t3perPos[i+j];
			avg_t5 += t5perPos[i+j];
		}
		avg_cons = avg_cons / window_size;
		avg_c3 /= window_size;
		avg_c5 /= window_size;
		avg_t3 /= window_size;
		avg_t5 /= window_size;

		InterestingRegion* temp;
		if (avg_c3 > 0){
			temp = new InterestingRegion(i, i+window_size-1,trueHelix, avg_c3, avg_cons);
			interestingRegions.push_back(temp);
		}
		if (avg_c5 > 0){
			temp = new InterestingRegion(i, i+window_size-1,trueHelix, avg_c5, avg_cons);
			interestingRegions.push_back(temp);
		}
		if (avg_t3 > 0){
			temp = new InterestingRegion(i, i+window_size-1,trueHelix, avg_t3, avg_cons);
			interestingRegions.push_back(temp);
		}
		if (avg_t5 > 0){
			temp = new InterestingRegion(i, i+window_size-1,trueHelix, avg_t5, avg_cons);
			interestingRegions.push_back(temp);
		}


	}



	//	interestingRegions[interestingRegions.size() - 1]->print();



}
string Alignment::competingHelicesHistogram(bool naive){
	vector<int> histo;
	stringstream ss;
	int pos;
	for (unsigned int j = 0; j < alignedSeqs[0]->size(); j++) {
		histo.push_back(0);
	}
	for (unsigned int seqIndex = 0; seqIndex < alignedSeqs.size(); seqIndex++){
		pos = 0;
		for (unsigned int j = 0; j < alignedSeqs[seqIndex]->size(); j++){
			if (alignedSeqs[seqIndex]->at(j) != '-' ) {
				int helixLabelatPos = helixLabels[seqIndex]->at(pos);
				for (int competingSeqIndex = 0; competingSeqIndex < (int)competingHelices[seqIndex]->size(); competingSeqIndex++) {
					//if (naive || (competingSeqIndex != helixLabelatPos) ) {
					//if (naive || (competingHelices[seqIndex]->at(competingSeqIndex)->trueHelixIndex != helixLabelatPos) ) {
					if (naive || StatsMatrix[seqIndex]->at(helixLabelatPos)->at(competingSeqIndex)->isZero() ) {
						//cout << "DEBUG " << helixLabels[i]->at(pos) << " , " << competingSeqIndex << endl;
						histo[j] += competingHelices[seqIndex]->at(competingSeqIndex)->isPosCompeting(pos);
						//						cout << "DEBUG-trueHelixIndex = " << competingHelices[seqIndex]->at(competingSeqIndex)->trueHelixIndex << endl;
					} else {
						//						cout << "DEBUG not adding " << competingHelices[seqIndex]->at(competingSeqIndex)->trueHelixIndex << " " << helixLabelatPos << endl;
					}
					//exit(0);
				}
				pos++;
			}
		}
	}
	for(unsigned int j = 0; j < histo.size(); j++) {
		ss << j << "\t" << histo[j] << "\n";
	}
	return ss.str();
}

string Alignment::competingHelicesHistogramTESTING(bool naive){
	vector<double> histo;
	stringstream ss;
	/* tmp rgoya */
	ofstream file_tmp("true_helix_histo/all-helices-ind-counts.txt");
	/* end tmp rgoya */
	for (unsigned int j = 0; j < alignedSeqs[0]->size(); j++) {
		histo.push_back(0);
	}
	vector<int> pos;
	vector<unsigned int> helixNumCompetingAtPos;
	for(unsigned int j = 0; j < numSeqs(); j++) {
		pos.push_back(0);
		helixNumCompetingAtPos.push_back(0);
	}

	for (unsigned int j = 0; j < alignedSeqs[0]->size(); j++){
		for (unsigned int seqIndex = 0; seqIndex < alignedSeqs.size(); seqIndex++){
			helixNumCompetingAtPos[seqIndex] = 0;
			if (alignedSeqs[seqIndex]->at(j) != '-' ) {
				int helixLabelatPos = helixLabels[seqIndex]->at(pos[seqIndex]);
				for (int competingSeqIndex = 0; competingSeqIndex < (int)competingHelices[seqIndex]->size(); competingSeqIndex++) {
					//if (naive || (competingHelices[seqIndex]->at(competingSeqIndex)->trueHelixIndex != helixLabelatPos) ) {
					if(helixLabelatPos == -1) {
						helixNumCompetingAtPos[seqIndex] += competingHelices[seqIndex]->at(competingSeqIndex)->isPosCompeting(pos[seqIndex]);
					} else if(naive || StatsMatrix[seqIndex]->at(helixLabelatPos)->at(competingSeqIndex)->isZero() ) {
						helixNumCompetingAtPos[seqIndex] += competingHelices[seqIndex]->at(competingSeqIndex)->isPosCompeting(pos[seqIndex]);
					}
					//exit(0);
				}
				pos[seqIndex]++;
			}
		}
		int tmpComp = 0;
		ss << j ;
		for(unsigned int i = 0; i < numSeqs(); i++) {
			if(helixNumCompetingAtPos[i] > 0)
				tmpComp++;
		}
		ss << "\t" <<  (double) tmpComp/numSeqs() ;
		/*
		ss << "DEBUG " << j ;
		for(unsigned int i = 0; i < numSeqs(); i++) {
			ss << "\t" << helixNumCompetingAtPos[i];
		}
		 */
		ss << "\n";

		/*
		 * tmp rgoya
		 * print output file with sequence counts
		 */
		file_tmp << j;
		for(unsigned int i = 0; i < numSeqs(); i++) {
			file_tmp << "\t" << helixNumCompetingAtPos[i];
		}
		file_tmp << "\n";
		/*
		 * end print
		 */
	}

	file_tmp.close();
	return ss.str();
}

/*
string Alignment::competingHelicesHistogramTESTING(bool naive) {
	stringstream ss;
	unsigned int totalSum = 0;

	vector<int> pos;
	vector<int> helixNumCompetingAtPos;
	vector<int> helixNumCompetingAtPos_c3;
	vector<int> helixNumCompetingAtPos_t3;
	vector<int> helixNumCompetingAtPos_c5;
	vector<int> helixNumCompetingAtPos_t5;
	for(unsigned int j = 0; j < numSeqs(); j++) {
		pos.push_back(0);
		helixNumCompetingAtPos.push_back(0);
		helixNumCompetingAtPos_c3.push_back(0);
		helixNumCompetingAtPos_t3.push_back(0);
		helixNumCompetingAtPos_c5.push_back(0);
		helixNumCompetingAtPos_t5.push_back(0);
	}

	for (unsigned int j = 0; j < alignedSeqs[0]->size(); j++){
		for (unsigned int seqIndex = 0; seqIndex < alignedSeqs.size(); seqIndex++){
			helixNumCompetingAtPos[seqIndex] = 0;
			helixNumCompetingAtPos_c3[seqIndex] = 0;
			helixNumCompetingAtPos_t3[seqIndex] = 0;
			helixNumCompetingAtPos_c5[seqIndex] = 0;
			helixNumCompetingAtPos_t5[seqIndex] = 0;
			if ( alignedSeqs[seqIndex]->at(j) != '-' ) {
				int helixLabelatPos = helixLabels[seqIndex]->at(pos[seqIndex]);
				for (int competingSeqIndex = 0; competingSeqIndex < competingHelices[seqIndex]->size(); competingSeqIndex++) {

					if( naive || StatsMatrix[seqIndex]->at(helixLabelatPos)->at(competingSeqIndex)->isZero()  ) {
						//if ( competingHelices[seqIndex]->at(competingSeqIndex)->trueHelixIndex == trueHelix ) {
						for(unsigned int trueHelix = 0; trueHelix < StatsMatrix[seqIndex]->size(); trueHelix++) {
							if (!StatsMatrix[seqIndex]->at(trueHelix)->at(competingSeqIndex)->isZero() ) {
								assert(helixNumCompetingAtPos[seqIndex] >= 0);
								if(competingHelices[seqIndex]->at(competingSeqIndex)->isPosCompeting(pos[seqIndex])) {
									//helixNumCompetingAtPos[seqIndex] += competingHelices[seqIndex]->at(competingSeqIndex)->isPosCompeting(pos[seqIndex]);
									helixNumCompetingAtPos[seqIndex]++;
									if(StatsMatrix[seqIndex]->at(trueHelix)->at(competingSeqIndex)->cis5() > 0) {
										//cout << "DEBUG cis5" << StatsMatrix[seqIndex]->at(trueHelix)->at(competingSeqIndex)->cis5() << endl;
										helixNumCompetingAtPos_c5[seqIndex]++;
									}
									if(StatsMatrix[seqIndex]->at(trueHelix)->at(competingSeqIndex)->cis3() > 0) {
										helixNumCompetingAtPos_c3[seqIndex]++;
									}
									if(StatsMatrix[seqIndex]->at(trueHelix)->at(competingSeqIndex)->trans5() > 0) {
										helixNumCompetingAtPos_t5[seqIndex]++;
									}
									if(StatsMatrix[seqIndex]->at(trueHelix)->at(competingSeqIndex)->trans3() > 0) {
										helixNumCompetingAtPos_t3[seqIndex]++;
									}
								}
							}
						}
					} else {
						assert(helixNumCompetingAtPos[seqIndex] == 0);
						helixNumCompetingAtPos[seqIndex] = -1;
					}
				}
				pos[seqIndex]++;
			}
		}
		int tmpTrue = 0;
		int tmpComp = 0;
		int c3 = 0;
		int c5 = 0;
		int t3 = 0;
		int t5 = 0;
		ss << j ;
		for(unsigned int i = 0; i < numSeqs(); i++) {
			if(helixNumCompetingAtPos[i] > 0)
				tmpComp++;
			else if(helixNumCompetingAtPos[i] < 0)
				tmpTrue--;
			if(helixNumCompetingAtPos_c3[i] > 0) { c3++; }
			if(helixNumCompetingAtPos_c5[i] > 0) { c5++; }
			if(helixNumCompetingAtPos_t3[i] > 0) { t3++; }
			if(helixNumCompetingAtPos_t5[i] > 0) { t5++; }
		}
		totalSum += tmpComp;
		ss << "\t" << (double) tmpTrue/numSeqs() << "\t" << (double) tmpComp/numSeqs() ;
		ss << "\t" << (double) c3/numSeqs();
		ss << "\t" << (double) c5/numSeqs();
		ss << "\t" << (double) t3/numSeqs();
		ss << "\t" << (double) t5/numSeqs();
		ss << "\n";
	}
	return ss.str();

}
 */
string Alignment::competingHelicesHistogramPerHelixTESTING(int trueHelix){
	stringstream ss;
	unsigned int totalSum = 0;

	vector<int> pos;
	vector<int> helixNumCompetingAtPos;
	vector<int> helixNumCompetingAtPos_c3;
	vector<int> helixNumCompetingAtPos_t3;
	vector<int> helixNumCompetingAtPos_c5;
	vector<int> helixNumCompetingAtPos_t5;
	for(unsigned int j = 0; j < numSeqs(); j++) {
		pos.push_back(0);
		helixNumCompetingAtPos.push_back(0);
		helixNumCompetingAtPos_c3.push_back(0);
		helixNumCompetingAtPos_t3.push_back(0);
		helixNumCompetingAtPos_c5.push_back(0);
		helixNumCompetingAtPos_t5.push_back(0);
	}
	for (unsigned int j = 0; j < alignedSeqs[0]->size(); j++){
		for (unsigned int seqIndex = 0; seqIndex < alignedSeqs.size(); seqIndex++){
			helixNumCompetingAtPos[seqIndex] = 0;
			helixNumCompetingAtPos_c3[seqIndex] = 0;
			helixNumCompetingAtPos_t3[seqIndex] = 0;
			helixNumCompetingAtPos_c5[seqIndex] = 0;
			helixNumCompetingAtPos_t5[seqIndex] = 0;
			if ( alignedSeqs[seqIndex]->at(j) != '-' ) {
				int helixLabelatPos = helixLabels[seqIndex]->at(pos[seqIndex]);
				if( trueHelix != helixLabelatPos ) {
					for (int competingSeqIndex = 0; competingSeqIndex < (int)competingHelices[seqIndex]->size(); competingSeqIndex++) {
						//if ( competingHelices[seqIndex]->at(competingSeqIndex)->trueHelixIndex == trueHelix ) {
						if (!StatsMatrix[seqIndex]->at(trueHelix)->at(competingSeqIndex)->isZero() ) {
							assert(helixNumCompetingAtPos[seqIndex] >= 0);
							if(competingHelices[seqIndex]->at(competingSeqIndex)->isPosCompeting(pos[seqIndex])) {
								//helixNumCompetingAtPos[seqIndex] += competingHelices[seqIndex]->at(competingSeqIndex)->isPosCompeting(pos[seqIndex]);
								helixNumCompetingAtPos[seqIndex]++;
								if(StatsMatrix[seqIndex]->at(trueHelix)->at(competingSeqIndex)->cis5() > 0) {
									//cout << "DEBUG cis5" << StatsMatrix[seqIndex]->at(trueHelix)->at(competingSeqIndex)->cis5() << endl;
									helixNumCompetingAtPos_c5[seqIndex]++;
								}
								if(StatsMatrix[seqIndex]->at(trueHelix)->at(competingSeqIndex)->cis3() > 0) {
									helixNumCompetingAtPos_c3[seqIndex]++;
								}
								if(StatsMatrix[seqIndex]->at(trueHelix)->at(competingSeqIndex)->trans5() > 0) {
									helixNumCompetingAtPos_t5[seqIndex]++;
								}
								if(StatsMatrix[seqIndex]->at(trueHelix)->at(competingSeqIndex)->trans3() > 0) {
									helixNumCompetingAtPos_t3[seqIndex]++;
								}
							}
						}
					}
				} else {
					assert(helixNumCompetingAtPos[seqIndex] == 0);
					helixNumCompetingAtPos[seqIndex] = -1;
				}
				pos[seqIndex]++;
			}
		}
		int tmpTrue = 0;
		int tmpComp = 0;
		int c3 = 0;
		int c5 = 0;
		int t3 = 0;
		int t5 = 0;
		for(unsigned int i = 0; i < numSeqs(); i++) {
			if(helixNumCompetingAtPos[i] > 0)
				tmpComp++;
			else if(helixNumCompetingAtPos[i] < 0)
				tmpTrue--;
			if(helixNumCompetingAtPos_c3[i] > 0) { c3++; }
			if(helixNumCompetingAtPos_c5[i] > 0) { c5++; }
			if(helixNumCompetingAtPos_t3[i] > 0) { t3++; }
			if(helixNumCompetingAtPos_t5[i] > 0) { t5++; }
		}
		totalSum += tmpComp;
		ss << j ;
		ss << "\t" << (double) tmpTrue/numSeqs() << "\t" << (double) tmpComp/numSeqs() ;
		ss << "\t" << (double) c3/numSeqs();
		ss << "\t" << (double) c5/numSeqs();
		ss << "\t" << (double) t3/numSeqs();
		ss << "\t" << (double) t5/numSeqs();
		ss << "\n";
	}

	cout << "True helix " << trueHelix << " with " << totalSum << " competing positions" << endl;
	return ss.str();
}
/*
vector<vector<vector<double> > > Alignment::competingHelicesHistogramPerHelixVector(int trueHelix){
	vector<vector<vector<double> > > histo;
	unsigned int totalSum = 0;

	vector<int> pos;
	vector<int> helixNumCompetingAtPos;
	vector<int> helixNumCompetingAtPos_c3;
	vector<int> helixNumCompetingAtPos_t3;
	vector<int> helixNumCompetingAtPos_c5;
	vector<int> helixNumCompetingAtPos_t5;
	for(unsigned int j = 0; j < numSeqs(); j++) {
		pos.push_back(0);
		helixNumCompetingAtPos.push_back(0);
		helixNumCompetingAtPos_c3.push_back(0);
		helixNumCompetingAtPos_t3.push_back(0);
		helixNumCompetingAtPos_c5.push_back(0);
		helixNumCompetingAtPos_t5.push_back(0);
	}
	for (unsigned int j = 0; j < alignedSeqs[0]->size(); j++){
		for (unsigned int seqIndex = 0; seqIndex < alignedSeqs.size(); seqIndex++){
			helixNumCompetingAtPos[seqIndex] = 0;
			helixNumCompetingAtPos_c3[seqIndex] = 0;
			helixNumCompetingAtPos_t3[seqIndex] = 0;
			helixNumCompetingAtPos_c5[seqIndex] = 0;
			helixNumCompetingAtPos_t5[seqIndex] = 0;
			if ( alignedSeqs[seqIndex]->at(j) != '-' ) {
				int helixLabelatPos = helixLabels[seqIndex]->at(pos[seqIndex]);
				if( trueHelix != helixLabelatPos ) {
					for (int competingSeqIndex = 0; competingSeqIndex < competingHelices[seqIndex]->size(); competingSeqIndex++) {
						//if ( competingHelices[seqIndex]->at(competingSeqIndex)->trueHelixIndex == trueHelix ) {
						if (!StatsMatrix[seqIndex]->at(trueHelix)->at(competingSeqIndex)->isZero() ) {
							assert(helixNumCompetingAtPos[seqIndex] >= 0);
							if(competingHelices[seqIndex]->at(competingSeqIndex)->isPosCompeting(pos[seqIndex])) {
							//helixNumCompetingAtPos[seqIndex] += competingHelices[seqIndex]->at(competingSeqIndex)->isPosCompeting(pos[seqIndex]);
								helixNumCompetingAtPos[seqIndex]++;
								if(StatsMatrix[seqIndex]->at(trueHelix)->at(competingSeqIndex)->cis5() > 0) {
									//cout << "DEBUG cis5" << StatsMatrix[seqIndex]->at(trueHelix)->at(competingSeqIndex)->cis5() << endl;
									helixNumCompetingAtPos_c5[seqIndex]++;
								}
								if(StatsMatrix[seqIndex]->at(trueHelix)->at(competingSeqIndex)->cis3() > 0) {
									helixNumCompetingAtPos_c3[seqIndex]++;
								}
								if(StatsMatrix[seqIndex]->at(trueHelix)->at(competingSeqIndex)->trans5() > 0) {
									helixNumCompetingAtPos_t5[seqIndex]++;
								}
								if(StatsMatrix[seqIndex]->at(trueHelix)->at(competingSeqIndex)->trans3() > 0) {
									helixNumCompetingAtPos_t3[seqIndex]++;
								}
							}
						}
					}
				} else {
					assert(helixNumCompetingAtPos[seqIndex] == 0);
					helixNumCompetingAtPos[seqIndex] = -1;
				}
				pos[seqIndex]++;
			}
		}
		int tmpTrue = 0;
		int tmpComp = 0;
		int c3 = 0;
		int c5 = 0;
		int t3 = 0;
		int t5 = 0;
		for(unsigned int i = 0; i < numSeqs(); i++) {
			if(helixNumCompetingAtPos[i] > 0)
				tmpComp++;
			else if(helixNumCompetingAtPos[i] < 0)
				tmpTrue--;
			if(helixNumCompetingAtPos_c3[i] > 0) { c3++; }
			if(helixNumCompetingAtPos_c5[i] > 0) { c5++; }
			if(helixNumCompetingAtPos_t3[i] > 0) { t3++; }
			if(helixNumCompetingAtPos_t5[i] > 0) { t5++; }
		}
		totalSum += tmpComp;
		histo.push_back();
		ss << j ;
		ss << "\t" << (double) tmpTrue/numSeqs() << "\t" << (double) tmpComp/numSeqs() ;
		ss << "\t" << (double) c3/numSeqs();
		ss << "\t" << (double) c5/numSeqs();
		ss << "\t" << (double) t3/numSeqs();
		ss << "\t" << (double) t5/numSeqs();
		ss << "\n";
	}

	cout << "True helix " << trueHelix << " with " << totalSum << " competing positions" << endl;
	return ss.str();
}
 */
unsigned int Alignment::numSeqs() {
	assert((unsigned int) alignedSeqs.size() == (unsigned int) competingHelices.size());
	assert((unsigned int) alignedSeqs.size() == (unsigned int) seqStructs.size());
	return (unsigned int) alignedSeqs.size();
}

unsigned int Alignment::numTrueHelices() {
	return helixNumber;
}

string Alignment::printCompetingHelices(int trueHelix) {
	stringstream ss;
	ss << "All Helices\n" << alignedStructString() << endl;
	ss << "True helix " << trueHelix << "\n";
	ss << printTrueHelix(trueHelix) << "\n";
	ss << "Competing helices\n";
	for(unsigned int seqIndex = 0; seqIndex < numSeqs(); seqIndex++) {
		ss << *seqNames[seqIndex] << "\n";
		for(unsigned int competingHelixIndex = 0; competingHelixIndex < StatsMatrix[seqIndex]->at(trueHelix)->size(); competingHelixIndex++) {
			if(StatsMatrix[seqIndex]->at(trueHelix)->at(competingHelixIndex)->trans3() > 0) {
				ss << competingHelices[seqIndex]->at(competingHelixIndex)->printAlignedHelix(*seq2AlignmentMap[seqIndex],alignedStruct.size());
				ss << "\t" << StatsMatrix[seqIndex]->at(trueHelix)->at(competingHelixIndex)->trans3() << "\n";
			}
		}
	}
	ss << "True helix " << trueHelix << "\n";
	ss << printTrueHelix(trueHelix) << "\n";
	return ss.str();
}

unsigned int Alignment::numCompetingOnTrueHelix(int trueHelix) {
	/* Using Stats matrix:
	 * 	[Sequence][True Helices][Competing Helices]
	 */
	unsigned int totalSum = 0;
	for(unsigned int seqIndex = 0; seqIndex < alignedSeqs.size(); seqIndex++) {
		for(unsigned int competingSeqIndex = 0; competingSeqIndex < StatsMatrix[seqIndex]->at(trueHelix)->size(); competingSeqIndex++) {
			if(!StatsMatrix[seqIndex]->at(trueHelix)->at(competingSeqIndex)->isZero()) {
				totalSum++;
			}
		}
	}
	return totalSum;
}

unsigned int Alignment::numCompetingOnSeqTrueHelix(int seqIndex, int trueHelix) {
	/* Using Stats matrix:
	 * 	[Sequence][True Helices][Competing Helices]
	 */
	unsigned int totalSum = 0;
	for(unsigned int competingSeqIndex = 0; competingSeqIndex < StatsMatrix[seqIndex]->at(trueHelix)->size(); competingSeqIndex++) {
		if(!StatsMatrix[seqIndex]->at(trueHelix)->at(competingSeqIndex)->isZero()) {
			totalSum++;
		}
	}
	return totalSum;
}



void Alignment::readAlignment(string & filename, bool includesStruct)
{
	ifstream fastaAln(filename.c_str());
	string line;
	string* seq = NULL;
	string* name = NULL;
	if(!fastaAln.is_open()){
		cerr << "Error: cannot open alignment file " << filename << endl;
		exit(-1);
	}
	while(fastaAln.good()){
		getline(fastaAln, line);
		Utilities::TrimSpaces(line);
		//skip blank lines
		if(line.length() == 0){
			continue;
		}
		if(line.at(0) == '>'){
			if(seq){
				//convert to upper case, and convert Ts to Us
				for(unsigned int i = 0; i < seq->length(); i++){
					seq->at(i) = toupper(seq->at(i));
					if(seq->at(i) == 'T'){
						seq->at(i) = 'U';
					}
				}
				alignedSeqs.push_back(seq);
			}
			line = line.substr(1);
			name = new string(line);
			seqNames.push_back(name);
			seq = new string();
		}
		else{

			*seq = *seq + line;
		}
	}

	assert(seq->length() > 0);
	if (includesStruct){
		//last sequence contains the structure
		seqNames.pop_back();
		parseDotBracket(*seq);
		delete seq;
	}
	else{
		//convert to upper case, and convert Ts to Us
		for(unsigned int i = 0; i < seq->length(); i++){
			seq->at(i) = toupper(seq->at(i));
			if(seq->at(i) == 'T'){
				seq->at(i) = 'U';
			}
		}
		alignedSeqs.push_back(seq);
	}

	fastaAln.close();

	//extract unaligned seqs and make map of unaligned seq positions to aligned seq positions
	vector<int>* temp;
	vector<int>* temp2;

	unsigned int alignedSeqLength = alignedSeqs[0]->length();

	for (unsigned int i = 0; i < alignedSeqs.size(); i++){

		assert(alignedSeqLength == alignedSeqs[i]->length()); //sanity check
		temp = new vector<int>();
		temp2 = new vector<int>();
		seq = new string();

		for (unsigned int j = 0; j < alignedSeqLength; j++){
			if (alignedSeqs[i]->at(j) != '-'){
				temp->push_back(j);
				*seq = *seq + alignedSeqs[i]->at(j);
				temp2->push_back(temp->size() - 1);
			}
			else{
				temp2->push_back(-1);
			}
		}
		seqs.push_back(seq);
		seq2AlignmentMap.push_back(temp);
		alignment2SeqMap.push_back(temp2);

	}

	if(includesStruct){
		fillSeqStructs();
	}

}

void Alignment::parseDotBracket(string & dotBracketString)
{
	//TODO: remove duplicate function in AlignmentGenerator class
	stack<int> parenStack, bracketStack, ltgtStack, braceStack, aStack, bStack, cStack, dStack;

	alignedStruct.clear();
	//sanity check
	assert(dotBracketString.length() == alignedSeqs[0]->length());
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

		alignedStruct.push_back(pairingPos);
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
}

void Alignment::readStruct(string & filename)
{


	ifstream ctFile(filename.c_str());
	string line;

	if(!ctFile.is_open()){
		cerr << "Error: cannot structure file " << filename << endl;
		exit(-1);
	}

	//skip first line:
	int spacerPos;
	getline(ctFile, line);

	while(ctFile.good()){
		getline(ctFile, line);
		Utilities::TrimSpaces(line);

		//skip blank lines
		if(line.length() == 0){
			continue;
		}
		//cout << line << endl;
		for(int i = 0; i < 4; i++){ //remove columns 1-4
			spacerPos = line.find_first_of("\t ");
			line = line.substr(spacerPos+1);
			Utilities::TrimSpaces(line); //removes leading spaces/tabs (if columns are separated by more than one).
		}
		//cout << line << endl;
		for(int i = 0; i < 1; i++){ //remove column 6
			spacerPos = line.find_last_of("\t ");
			line = line.substr(0, spacerPos);
			Utilities::TrimSpaces(line); //removes trailing spaces/tabs (if columns are separated by more than one).
		}
		//cout << line << endl;
		istringstream buffer(line);
		int temp;
		buffer >> temp;

		alignedStruct.push_back(temp -1);
	}

	fillSeqStructs();

}

void Alignment::fillSeqStructs(){

	assert(alignedStruct.size() == alignedSeqs[0]->length());

	//clear seqStructs:

	for(unsigned int i = 0; i < seqStructs.size(); i++){
			delete seqStructs[i];
	}
	seqStructs.clear();


	//fill vectors for structure per sequence
	vector<int>* tempVec;
	char a, b;

	for(unsigned int i = 0; i < alignedSeqs.size();i++){
		tempVec = new vector<int>();
		for(unsigned int j = 0; j < alignedSeqs[i]->length(); j++){
			a = alignedSeqs[i]->at(j);
			if (!(a == '-')){
				if (alignedStruct[j] == -1){
					tempVec->push_back(alignedStruct[j]); //push back -1
				}
				else{
					b = alignedSeqs[i]->at(alignedStruct[j]);
					if(Utilities::validBP(a,b)){
						tempVec->push_back(alignment2SeqMap[i]->at(alignedStruct[j]));
					}
					else{
						tempVec->push_back(-1);
					}
				}
			}
			//else do nothing -> we are at a gap
		}
		assert(tempVec->size() == seqs[i]->length());
		seqStructs.push_back(tempVec);

	}

}

void Alignment::alignmentGStats(const double eStack[4][4][4][4])
{
	for(unsigned int i = 0; i < seqs.size(); i++){
		sequenceGStats(eStack, i);
	}

}

void Alignment::sequenceGStats(const double eStack[4][4][4][4], const int seqIndex)
{
	int **DP;
	int i, j, k, i_2, j_2;
	double eE;

	string Seq = *(seqs[seqIndex]);
	int iL = Seq.length();

	vector<int> P = *(seqStructs[seqIndex]);
	// Initialize values

	//info needed to store competing helices
	int start5, end5, start3, end3;
	set<int> competedHelices;
	set<int>::iterator cursor;
	CompetingHelix* currentCompetingHelix;
	int helixIndex;

	//clear StatsMatrix for current sequence
	for(unsigned int index = 0; index < StatsMatrix[seqIndex]->size(); index++){
		for(unsigned int index2 = 0; index2 < StatsMatrix[seqIndex]->at(index)->size(); index2++){
			delete StatsMatrix[seqIndex]->at(index)->at(index2);
		}
		StatsMatrix[seqIndex]->at(index)->clear();
	}
	for(unsigned int index = 0; index < competingHelices[seqIndex]->size(); index++){
		delete competingHelices[seqIndex]->at(index);
	}
	competingHelices[seqIndex]->clear();

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

			if ((DP[i][j] == 0) && (DP[i + 1][j - 1] > minStemLength)) {
				// iterating on one found helix


				if(isHelixCompeting(seqIndex, i+1, j-1, DP[i+1][j-1])){
					//have a competing helix!
					//cout << "competing helix found\n";
					start5 = i + 1;
					end3 = j - 1;

					//add a row of StatsWrappers to the StatsMatrix for sequence seqIndex
					for(unsigned int index = 0; index < StatsMatrix[seqIndex]->size(); index++){
						StatsMatrix[seqIndex]->at(index)->push_back(new StatsWrapper());
					}

					eE = 0.0;

					for (k = 0; k < DP[i + 1][j - 1] - 1; k++) {

						i_2 = i + 1 + k;
						j_2 = j - 1 - k;

						eE += eStack[Utilities::iAt(Seq[i_2])][Utilities::iAt(Seq[i_2 + 1])][Utilities::iAt(Seq[j_2 - 1])][Utilities::iAt(Seq[j_2])];
					}
					//cout << "eE: " << eE << endl;
					//instead, calculate p stats ----
					//eE = 1;
					//-------------------------------
					for (k = 0; k < DP[i + 1][j - 1]; k++) {

						i_2 = i + 1 + k;
						j_2 = j - 1 - k;

						// if i paired

						if ((P[i_2] != -1) && (Utilities::iCheckLongHelix(P, i_2) == 1)) {


							helixIndex = helixLabels[seqIndex]->at(i_2);

							if (P[i_2] < i_2) { // (p,i,j) --> (i^, i, c) 3'cis

								double statvalue = -(eE) / ((double) (j - 2 - 2 * k - i) * log((double) (iL - i_2)));

								StatsMatrix[seqIndex]->at(helixIndex)->back()->eStats[0] += statvalue;
								StatsMatrix[seqIndex]->at(helixIndex)->back()->iStats[0]++;

							}
							else if ((P[i_2] > i_2) && (P[i_2] < j_2)) { // (i,p,j) --> (i, i^, c) 3'trans

								double statvalue = -(eE) / ((double) (j_2 - P[i_2]) * log((double) (iL - P[i_2])));

								StatsMatrix[seqIndex]->at(helixIndex)->back()->eStats[1] += statvalue;
								StatsMatrix[seqIndex]->at(helixIndex)->back()->iStats[1]++;

							}
							else if (P[i_2] > j_2) { // (i,j,p) --> (i, c, i^) 3' mid

								double statvalue = -(eE) / ((double) (j - 2 - 2 * k - i)
										* log((double) (P[i_2] - i_2)));
								StatsMatrix[seqIndex]->at(helixIndex)->back()->eStats[2] += statvalue;
								StatsMatrix[seqIndex]->at(helixIndex)->back()->iStats[2]++;
							}
						}

						// if j paired

						if ((P[j_2] != -1) && (Utilities::iCheckLongHelix(P, j_2) == 1)) {

							//get helix label of true helix that is being competed
							helixIndex = helixLabels[seqIndex]->at(j_2);

							if (P[j_2] > j_2) { //(i,j,p) --> (c, i, i^) 5'cis

								double statvalue = -(eE) / ((double) (j - 2 - 2 * k
										- i) * log((double) (j_2)));

								StatsMatrix[seqIndex]->at(helixIndex)->back()->revEStats[0] += statvalue;
								StatsMatrix[seqIndex]->at(helixIndex)->back()->revIStats[0]++;

							}
							else if ((P[j_2] > i_2) && (P[j_2] < j_2)) { //(i,p,j) --> (c, i^, i) 5'trans

								double statvalue = -(eE) / ((double) (P[j_2] - i_2) * log((double) (P[j_2])));
								StatsMatrix[seqIndex]->at(helixIndex)->back()->revEStats[1] += statvalue;
								StatsMatrix[seqIndex]->at(helixIndex)->back()->revIStats[1]++;
							}
							else if (P[j_2] < i_2) { //(p, i, j) --> (i^, c, i) 5' mid
								double statvalue = -(eE) / ((double) (i_2 - P[j_2]) * log((double) (j_2 - P[j_2])));
								StatsMatrix[seqIndex]->at(helixIndex)->back()->revEStats[2] += statvalue;
								StatsMatrix[seqIndex]->at(helixIndex)->back()->revIStats[2]++;
							}
						}
					}

					end5 = i + 1 + k - 1;
					start3 = j - 1 - k + 1;

					currentCompetingHelix = new CompetingHelix(start5, end5, start3, end3);
					competingHelices[seqIndex]->push_back(currentCompetingHelix);

				}//if (isHelixCompeting)
				//if helix is not competing, ignore it
			}
		}
	}


	for (i = 0; i < iL; i++) {
		free(DP[i]);
	}
	free(DP);

	return;
}

//bool Alignment::validBP(char a, char b)
//{
//
//	b = toupper(b);
//	a = toupper(a);
//
//	//note: I think a lot of other methods are not safe for Ts
//	if (a == 'T'){
//		a = 'U';
//	}
//	else if(b == 'T'){
//		b = 'U';
//	}
//
//	string temp("");
//	temp += a;
//	temp += b;
//
//	if (temp.compare("AU") == 0 || temp.compare("UA") == 0){
//		return true;
//	}
//
//	if (temp.compare("GU") == 0 || temp.compare("UG") == 0){
//		return true;
//	}
//
//	if (temp.compare("GC") == 0 || temp.compare("CG") == 0){
//		return true;
//	}
//	//cout << "failed: " << temp << endl;
//	return false;
//}



void Alignment::labelHelices2()
{
	//clear data
	alignmentLabels.clear();
	for(unsigned int i = 0; i < helixLabels.size(); i++){
		delete helixLabels[i];
	}
	helixLabels.clear();

	for(unsigned int i = 0; i < trueHelices.size(); i++){
		delete trueHelices[i];
	}
	trueHelices.clear(); //clears trueHelices, but doesn't refill them

	//set up empty vector
	for (unsigned int i = 0; i < alignedStruct.size(); i++) {
		alignmentLabels.push_back(-1);
	}

	for (unsigned int j = 0; j < seqStructs.size(); j++){
		vector<int>* temp = new vector<int>();
		for (unsigned int i = 0; i < seqStructs[j]->size(); i++){
			temp->push_back(-1);
		}
		helixLabels.push_back(temp);
	}

	int helixCounter = -1;
	//populate labels
	for (int i = 0; i < (int)alignedStruct.size(); i++) {
		if (alignedStruct[i] == -1) {
			//do nothing
		} else if (alignedStruct[i] > i) {
			//start of a helix:
			helixCounter++;
			alignmentLabels.at(i) = helixCounter;
			alignmentLabels.at(alignedStruct[i]) = helixCounter;

			//set labels for each sequence
			for(unsigned int j = 0; j < helixLabels.size(); j++){
				if(alignment2SeqMap[j]->at(i) != -1 && alignment2SeqMap[j]->at(alignedStruct[i]) != -1
						&& seqStructs[j]->at(alignment2SeqMap[j]->at(i)) != -1
						&& seqStructs[j]->at(alignment2SeqMap[j]->at(alignedStruct[i])) != -1){
					helixLabels[j]->at(alignment2SeqMap[j]->at(i)) = helixCounter;
					helixLabels[j]->at(alignment2SeqMap[j]->at(alignedStruct[i])) = helixCounter;
				}
			}

			//extend helix as far as possible:
			i++;
			while(true){

				while (alignedStruct[i] == alignedStruct[i - 1] - 1) {
					alignmentLabels.at(i) = helixCounter;
					alignmentLabels.at(alignedStruct[i]) = helixCounter;

					for(unsigned int j = 0; j < helixLabels.size(); j++){
						if(alignment2SeqMap[j]->at(i) != -1 && alignment2SeqMap[j]->at(alignedStruct[i]) != -1
								&& seqStructs[j]->at(alignment2SeqMap[j]->at(i)) != -1
								&& seqStructs[j]->at(alignment2SeqMap[j]->at(alignedStruct[i])) != -1){
							//checks that, in this sequence (j), neither the paired bases are a gap and that they are paired
							helixLabels[j]->at(alignment2SeqMap[j]->at(i)) = helixCounter;
							helixLabels[j]->at(alignment2SeqMap[j]->at(alignedStruct[i])) = helixCounter;
						}
					}
					i++;
				}

				int endCursor = alignedStruct[i-1] - 1;

				//find the next bp:
				while(alignedStruct[i] == -1){
					i++;
				}

				if(alignedStruct[i] > i){
					while(alignedStruct[endCursor] == -1){
						endCursor--;
					}
					if(alignedStruct[i] == endCursor){
						//helices are continuous - keep labelling them the same

						alignmentLabels.at(i) = helixCounter;
						alignmentLabels.at(alignedStruct[i]) = helixCounter;

						//set labels for each sequence
						for(unsigned int j = 0; j < helixLabels.size(); j++){
							if(alignment2SeqMap[j]->at(i) != -1 && alignment2SeqMap[j]->at(alignedStruct[i]) != -1
									&& seqStructs[j]->at(alignment2SeqMap[j]->at(i)) != -1
									&& seqStructs[j]->at(alignment2SeqMap[j]->at(alignedStruct[i])) != -1){
								helixLabels[j]->at(alignment2SeqMap[j]->at(i)) = helixCounter;
								helixLabels[j]->at(alignment2SeqMap[j]->at(alignedStruct[i])) = helixCounter;
							}
						}
						i++;
						continue;
					}
				}
				break;


			} //end of while(true)
			i--;
		}
	}

	helixNumber = helixCounter + 1;

}


/**
 * Labels each helix with an integer number >0, and prints out a string representation of the
 * labeled sequence
 */
void Alignment::labelHelices()
{

	//clear data
	alignmentLabels.clear();
	for(unsigned int i = 0; i < helixLabels.size(); i++){
		delete helixLabels[i];
	}
	helixLabels.clear();

	for(unsigned int i = 0; i < trueHelices.size(); i++){
		delete trueHelices[i];
	}
	trueHelices.clear();

	string helixString = string("");
	string bracketString = string("");

	//set up empty vector
	for (unsigned int i = 0; i < alignedStruct.size(); i++) {
		alignmentLabels.push_back(-1);
	}

	for (unsigned int j = 0; j < seqStructs.size(); j++){
		vector<int>* temp = new vector<int>();
		for (unsigned int i = 0; i < seqStructs[j]->size(); i++){
			temp->push_back(-1);
		}
		helixLabels.push_back(temp);
	}

	int helixCounter = -1;
	//populate labels
	for (int i = 0; i < (int)alignedStruct.size(); i++) {
		if (alignedStruct[i] == -1) {
			//do nothing
		} else if (alignedStruct[i] > i) {
			//start of a helix:
			helixCounter++;
			alignmentLabels.at(i) = helixCounter;
			alignmentLabels.at(alignedStruct[i]) = helixCounter;

			//set labels for each sequence
			for(unsigned int j = 0; j < helixLabels.size(); j++){
				if(alignment2SeqMap[j]->at(i) != -1 && alignment2SeqMap[j]->at(alignedStruct[i]) != -1
						&& seqStructs[j]->at(alignment2SeqMap[j]->at(i)) != -1
						&& seqStructs[j]->at(alignment2SeqMap[j]->at(alignedStruct[i])) != -1){
					helixLabels[j]->at(alignment2SeqMap[j]->at(i)) = helixCounter;
					helixLabels[j]->at(alignment2SeqMap[j]->at(alignedStruct[i])) = helixCounter;
				}
			}

			//create new true helix:
			Helix * trueHelix= new Helix(i, alignedStruct[i], 1);

			//extend helix as far as possible:
			i++;
			//unsigned int skip_j = 1;
			//unsigned int skip_i = 1;
			//while (alignedStruct[i] == alignedStruct[i - skip_i] - skip_j) {
			while (alignedStruct[i] == alignedStruct[i - 1] - 1) {
				alignmentLabels.at(i) = helixCounter;
				alignmentLabels.at(alignedStruct[i]) = helixCounter;

				for(unsigned int j = 0; j < helixLabels.size(); j++){
					if(alignment2SeqMap[j]->at(i) != -1 && alignment2SeqMap[j]->at(alignedStruct[i]) != -1
							&& seqStructs[j]->at(alignment2SeqMap[j]->at(i)) != -1
							&& seqStructs[j]->at(alignment2SeqMap[j]->at(alignedStruct[i])) != -1){
						//checks that, in this sequence (j), neither the paired bases are a gap and that they are paired
						helixLabels[j]->at(alignment2SeqMap[j]->at(i)) = helixCounter;
						helixLabels[j]->at(alignment2SeqMap[j]->at(alignedStruct[i])) = helixCounter;
					}
				}

				trueHelix->length++;
				i++;
				//check if there's a bulge
				//if(alignedStruct[i+1] == alignedStruct[i - skip_i] - skip_j) {  i++; skip_i++; }
				//if(alignedStruct[i] == alignedStruct[i - skip_i] - skip_j - 1) {  skip_j++; }

			}
			i--;
			trueHelices.push_back(trueHelix); //trueHelices[helixCounter] should equal trueHelix*

		}
	}

	helixNumber = helixCounter + 1;

	//Print Bracket Notation:
	/*
	for (int i = 0; i < alignedStruct.size(); i++) {
		if (P[i] > i) {
			bracketString += "(";
		} else if (P[i] == -1) {
			bracketString += ".";
		} else if (P[i] < i) {
			bracketString += ")";
		}
	}
	cout << bracketString << endl;
	 */
	//Print Labels
	/*
	for (int i = 0; i < iL; i++) {

		cout << PLabels.at(i);
	}
	cout << endl;

	cout << helixCounter << " helices found.\n";
	 */
}

void Alignment::printCompetingHelixCount()
{
	for(unsigned int i = 0; i < competingHelices.size(); i++){
		cout << competingHelices[i]->size() << endl;
	}
}

void Alignment::sequenceConservation()
{
	double baseCounts[4];
	char base;
	unsigned int seqIndex;
	double max;

	assert(!alignedSeqs.empty());
	for(unsigned int i = 0; i < alignedSeqs[0]->length(); i++){
		//loop over all positions in the alignment

		//set counts to zero
		for(int j = 0; j < 4; j++){
			baseCounts[j] = 0.0;
		}

		for(seqIndex = 0; seqIndex < alignedSeqs.size(); seqIndex++){
			base = alignedSeqs[seqIndex]->at(i);
			if (base == '-' || base == 'N'){
				//if position is gap, do nothing
			}
			else{
				baseCounts[Utilities::iAt(base)] += 1;
			}

		}
		//divide by the number of sequences to get
		//probability for each base
		for(int j = 0; j < 4; j++){
			baseCounts[j] = baseCounts[j] / alignedSeqs.size();
		}

		max = 0;
		for (int j = 0; j < 4; j++){
			if (baseCounts[j] > max){
				max = baseCounts[j];
			}
		}
		cout << max << endl;

	}

}

void Alignment::printCisTransGValuesPerHelix()
{
	//print header line
	for (int i = 0; i < helixNumber; i++){
		cout << "cis " << (i+1) << "\t" << "trans " << (i+1) << "\t";
	}
	cout << endl;


	for(unsigned int seqIndex = 0; seqIndex < seqs.size(); seqIndex++){
		for(int i = 0; i < helixNumber; i++){
			StatsWrapper sum = StatsWrapper(); //empty StatsWrapper
			for(unsigned int j = 0; j < StatsMatrix[seqIndex]->at(i)->size(); j++){
				sum += *(StatsMatrix[seqIndex]->at(i)->at(j));
			}
			cout << sum.cis() << "\t" << sum.trans() << "\t";
		}
		cout << "\n";
	}
}

void Alignment::printCisTransGValuesPerSequence()
{


	for(unsigned int seqIndex = 0; seqIndex < seqs.size(); seqIndex++){
		cout << "For comparison against get_all_statistics\n";
		cout << "Sequence " << seqIndex << " (" << *seqNames[seqIndex] << ")";
		StatsWrapper sum = StatsWrapper(); //empty StatsWrapper
		for(int i = 0; i < helixNumber; i++){
			for(unsigned int j = 0; j < StatsMatrix[seqIndex]->at(i)->size(); j++){
				sum += *(StatsMatrix[seqIndex]->at(i)->at(j));
			}
		}
		cout << "\t3cis= " << sum.cis3();
		cout << "\t5cis= " << sum.cis5();
		cout << "\t3trans= " << sum.trans3();
		cout << "\t5trans= " << sum.trans5();
		//cout << "\tcis= " << sum.cis() << "\ttrans= " << sum.trans() << "\n";
		cout << "\n";
	}
}

string Alignment::printTrueHelix(int helixIndex)
{
	string out = "";
	int pos = 0;
	for(unsigned int i = 0; i < alignmentLabels.size(); i++){
		if (helixIndex == alignmentLabels[i]){
			if (alignedStruct[i] > pos){
				out += "(";
			}
			else{
				out += ")";
			}
		}
		else{
			out += ".";
		}

		pos++;
	}

	return out;
}

void Alignment::printCisTransValues()
{
	//print header line
	cout << "cis\ttrans\n";

	for(unsigned int seqIndex = 0; seqIndex < seqs.size(); seqIndex++){
		StatsWrapper sum = StatsWrapper(); //empty StatsWrapper
		//loop over true helices
		for(int i = 0; i < helixNumber; i++){
			//loop over competing helices
			for(unsigned int j = 0; j < StatsMatrix[seqIndex]->at(i)->size(); j++){
				sum += *(StatsMatrix[seqIndex]->at(i)->at(j));
			}
		}
		cout << sum.cis() << "\t" << sum.trans() <<"\n";
	}
}



bool Alignment::isHelixCompeting(int seqIndex, int start_i, int start_j, int length)
{
	int i, j;
	bool competing = false;
	for (int k = 0; k < length; k++) {
		i = start_i + k;
		j = start_j - k;


		if(seqStructs[seqIndex]->at(i) == j){
			//i-j paired in the true structure
			//this helix is part of the true structure (or at least partially part of the true sequence)
			return false;
		}
		else if(seqStructs[seqIndex]->at(i) != -1){
			//helix competes with at least one true helix!
			competing = true;
		}
		else if(seqStructs[seqIndex]->at(j) != -1){
			//helix competes with at least one true helix!
			competing = true;
		}
	}
	//cout << "not competing\n";
	return competing;
}

void Alignment::alignmentInformationContent(){

	double baseCounts[4];
	char base;
	unsigned int seqIndex;
	double info;

	assert(!alignedSeqs.empty());

	for(unsigned int i = 0; i < alignedSeqs[0]->length(); i++){
		//loop over all positions in the alignment

		//set counts to zero
		for(int j = 0; j < 4; j++){
			baseCounts[j] = 0.0;
		}

		for(seqIndex = 0; seqIndex < alignedSeqs.size(); seqIndex++){
			base = alignedSeqs[seqIndex]->at(i);
			if (base == '-'){
				//if position is gap, treat as equal probability of being
				//any of the 4 nucleotides
				//this may not be the best way to go about this, but it
				//seems vaguely reasonable
				for(int j = 0; j < 4; j++){
					baseCounts[j] += 0.25;
				}
			}
			else{
				baseCounts[Utilities::iAt(base)] += 1;
			}

		}
		//divide by the number of sequences to get
		//probability for each base
		for(int j = 0; j < 4; j++){
			baseCounts[j] = baseCounts[j] / seqIndex; //after for loop, seqIndex == alignedSeqs.size()
		}

		//information content = log(4) - (- sum(P(base) * log(P(base)) ) )
		info = 2;
		for(int j = 0; j < 4; j++){
			if(baseCounts[j] > 0){
				info += baseCounts[j] * log(baseCounts[j]);
			} //important because log(0) is undefined -> for our purposes, 0 * log(0) = 0
		}
		//note:
		//http://helix.mcmaster.ca/721/outline2/node58.html and
		//Schneider et al. (1986) 'Information content of binding sites'
		//suggest using a correction factor to take into account
		//the finite sample size used to estimate P(base). I couldn't find
		//a straighforward description of how to calculate this correction
		//factor, however.
		//cout << (i+1) << "\t" << info << endl;
		cout << i << "\t" << info << endl;

	}
}

void  Alignment::printCompetingHelixMidpoints(){

	vector<int> seqLenths;
	for (unsigned int i = 0; i < seqs.size(); i++){
		seqLenths.push_back(seqs[i]->length());
	}

	for (unsigned int i = 0; i < competingHelices.size(); i++){
		for(unsigned int j = 0; j < competingHelices[i]->size(); j++){
			double mid = (competingHelices[i]->at(j)->end3 + competingHelices[i]->at(j)->start5) / 2.0;

			mid = mid / seqLenths[i];
			cout << mid << endl;
		}
	}

}

void Alignment::printCoreLogLikelihoods(Tree & treeRoot){
	if(printHeaders){
			cout << "fraction index\tfraction\tlength\tlogRatio\t3cis\t3trans\t5cis\t5trans\n";
		}
	//loop through core cutoffs:
	for(unsigned int i = 0; i < cores.size(); i++){

		if(cores[i]->size() > 0){
			double seqFrac = i; //core is present in this fraction of sequences of the alignment
			seqFrac /= seqNames.size();
			//cout << "core cutoff: " << i <<"\tFraction of seqs: " << seqFrac << endl;

			//loop through cores
			for (unsigned int j = 0; j < cores[i]->size();j++){

				double pairedLogLike = 0;
				double singleLogLike = 0;
				int is3cis=0;
				int is5cis=0;
				int is3trans=0;
				int is5trans=0;
				int pos5 = cores[i]->at(j)->getPos5();
				int pos3 = cores[i]->at(j)->getPos3();
				for (int k = 0; k < cores[i]->at(j)->getLength(); k++){
					pairedLogLike += log(treeRoot.calcFelsDouble(*this, pos5+k, pos3-k))/log(2);
					singleLogLike += log(treeRoot.calcFelsSingle(*this, pos5+k))/log(2);
					singleLogLike += log(treeRoot.calcFelsSingle(*this, pos3-k))/log(2);
					if(alignedStruct[pos5+k] != -1){
						if(alignedStruct[pos5+k] > pos5+k){
							is3trans++;
						}
						else{
							is3cis++;
						}
					}
					if(alignedStruct[pos3-k] != -1){
						if(alignedStruct[pos3-k] > pos3-k){
							is5cis++;
						}
						else{
							is5trans++;
						}
					}
				}




				double logRatio = (pairedLogLike - singleLogLike) / cores[i]->at(j)->getLength();
				cout << i << "\t" << seqFrac << "\t" << cores[i]->at(j)->getLength() <<"\t"<< logRatio << "\t"
				<< is3cis <<"\t" << is3trans <<"\t" << is5cis << "\t" << is5trans << endl;


				//cout << cores[i]->at(j)->dotBracket(alignedStruct.size()) << endl;
				//cout << cores[i]->at(j)->printAlignmentAtCore(alignedSeqs);
			}
		}
	}
}

void Alignment::printTrueHelixLogLikelihoods(Tree & treeRoot){
	if(printHeaders){
		cout << "Helix Index\tLength\tMidpoint\tlogRatio\n";
	}
	for(unsigned int i = 0; i < trueHelices.size(); i++){
		double pairedLogLike = 0;
		double singleLogLike = 0;
		int pos5 = trueHelices[i]->pos5;
		int pos3 = trueHelices[i]->pos3;
		for (int k = 0; k < trueHelices[i]->length; k++){
			pairedLogLike += log(treeRoot.calcFelsDouble(*this, pos5+k, pos3-k))/log(2);
			singleLogLike += log(treeRoot.calcFelsSingle(*this, pos5+k))/log(2);
			singleLogLike += log(treeRoot.calcFelsSingle(*this, pos3-k))/log(2);
		}

		double logRatio = (pairedLogLike - singleLogLike) / trueHelices[i]->length;
		cout << i << "\t" << trueHelices[i]->length <<"\t"
		<< trueHelices[i]->midpoint(alignedSeqs[0]->length())<< "\t" << logRatio <<  endl;

	}
}

void Alignment::printCompetingHelixLogLikelihoods(Tree & treeRoot){
	if(printHeaders){
		cout << "Sequence Index\tLength\tMidpoint\tlogRatio\tcorrectedLogRatio\tcis\ttrans\tmid\n";
	}
	for(unsigned int i = 0; i < competingHelices.size(); i++){
		//loop through competing helices on sequence i
		for(unsigned int j = 0; j < competingHelices[i]->size(); j++){
			//check if competing helix is 'clean'
//			int cleanCount = 0;
//			for(unsigned int k = 0; k < StatsMatrix[i]->size(); k++){
//				if(!(StatsMatrix[i]->at(k)->at(j)->isZero())){
//					cleanCount++;
//				}
//			}
//			if(cleanCount == 0){
//				cout << structString();
//				cout << competingHelices[i]->at(j)->printHelix2(seqs[i]->length()) << endl;
//
//			}
//			assert(cleanCount > 0);
			if(true){
			//only print 'clean' helices
			//if(helixIsClean(*(competingHelices[i]->at(j)), i)){

				double pairedLogLike = 0;
				double singleLogLike = 0;
				int pos5 = competingHelices[i]->at(j)->start5;
				int pos3 = competingHelices[i]->at(j)->end3;

				int length = competingHelices[i]->at(j)->end5 - pos5 + 1;
				//convert to alignment positions:
				for (int k = 0; k < length; k++){
					int alignedPos5 = seq2AlignmentMap[i]->at(pos5+k);
					int alignedPos3 = seq2AlignmentMap[i]->at(pos3-k);
					pairedLogLike += log(treeRoot.calcFelsDouble(*this, alignedPos5 , alignedPos3))/log(2);
					singleLogLike += log(treeRoot.calcFelsSingle(*this, alignedPos5))/log(2);
					singleLogLike += log(treeRoot.calcFelsSingle(*this, alignedPos3))/log(2);
				}

				double logRatio = pairedLogLike - singleLogLike;

				double midpoint = competingHelices[i]->at(j)->midpoint();
				midpoint = midpoint / seqs[i]->length();

				double cis = 0;
				double trans = 0;
				double mid = 0;
				//loop through true helices in stats matrix to find total cis and trans values
				for(unsigned int k = 0; k < StatsMatrix[i]->size();k++){
					cis += StatsMatrix[i]->at(k)->at(j)->cis();
					trans += StatsMatrix[i]->at(k)->at(j)->trans();
					mid += StatsMatrix[i]->at(k)->at(j)->mid();
				}
				//if(logRatio / (double) length > 0){
					cout << i << "\t" << length << "\t" << midpoint <<"\t"<< logRatio << "\t" << logRatio / (double)length
					<< "\t" << cis << "\t" << trans << "\t" << mid << endl;

					//cout << competingHelices[i]->at(j)->printHelix2(seqStructs[i]->size()) << endl;
					//cout << competingHelices[i]->at(j)->printAlignedHelix(*seq2AlignmentMap[i], alignedStruct.size()) << endl;
				//}
			}
		}
	}
}

void Alignment::printCompetingHelixTrueHelixMidpoints(Tree & treeRoot, const double eStack[4][4][4][4]){
	//print column names:
	cout << "Seq Index\tTrue Helix Index\tCompeting Helix Index\t"
	<<"True Helix Midpoint\tCompeting Helix Midpoint\t"
	<< "True Helix Length\tTrue Helix Log Ratio\t"
	<< "CompetingHelix Length\tCompeting Helix Log Ratio\t"
	<< "5'cis\t5'trans\t5'mid\t3'cis\t3'trans\t3'mid\t"
	<< "True Helix G\tCompetingHelixG"
	<< endl;
	//loop through sequences
	// i = sequence index
	for(unsigned int i = 0; i < StatsMatrix.size(); i++){
		//loop through true helices
		//j = true helix index
		for (unsigned int j = 0; j < StatsMatrix[i]->size(); j++){

			//figure out log likelihood for true helix:
			double pairedLogLike = 0;
			double singleLogLike = 0;
			int pos5 = trueHelices[j]->pos5;
			int pos3 = trueHelices[j]->pos3;
			for (int k = 0; k < trueHelices[j]->length; k++){
				pairedLogLike += log(treeRoot.calcFelsDouble(*this, pos5+k, pos3-k));
				singleLogLike += log(treeRoot.calcFelsSingle(*this, pos5+k));
				singleLogLike += log(treeRoot.calcFelsSingle(*this, pos3-k));
			}

			double trueHelixLogRatio = (pairedLogLike - singleLogLike) / trueHelices[j]->length;
			double trueHelixMidpoint = trueHelices[j]->midpoint(alignedSeqs[i]->length());

			//loop through competing helices
			//k = competing helix index
			for(unsigned int k = 0; k < StatsMatrix[i]->at(j)->size(); k++){
				StatsWrapper temp = *(StatsMatrix[i]->at(j)->at(k));
				//ignore this competing helix/true helix pair if they don't actually interact
				if(temp.isZero()){
					continue;
				}
				pairedLogLike = 0;
				singleLogLike = 0;

				pos5 = competingHelices[i]->at(k)->start5;
				pos3 = competingHelices[i]->at(k)->end3;

				//use alignment positions to map midpoint -> more straightforward than going
				//from alignment positions to sequence positions, since a sequence position
				//always has a position in the alignment, whereas an alignment position
				//might be a gap in a specific sequence.
				double competingHelixMidpoint = (seq2AlignmentMap[i]->at(pos5) + seq2AlignmentMap[i]->at(pos3)) / (2.0 * alignedSeqs[i]->length());

				int length = competingHelices[i]->at(k)->end5 - pos5 + 1;
				//convert to alignment positions:
				for (int l = 0; l < length; l++){
					int alignedPos5 = seq2AlignmentMap[i]->at(pos5+l);
					int alignedPos3 = seq2AlignmentMap[i]->at(pos3-l);
					pairedLogLike += log(treeRoot.calcFelsDouble(*this, alignedPos5 , alignedPos3));
					singleLogLike += log(treeRoot.calcFelsSingle(*this, alignedPos5));
					singleLogLike += log(treeRoot.calcFelsSingle(*this, alignedPos3));
				}

				double competingHelixLogRatio = (pairedLogLike - singleLogLike) / length;

				cout << i << "\t" << j << "\t" << k << "\t"
				<< trueHelixMidpoint << "\t" << competingHelixMidpoint << "\t"
				<< trueHelices[j]->length << "\t" << trueHelixLogRatio <<"\t"
				<< length << "\t" << competingHelixLogRatio << "\t"
				<< temp.cis5() << "\t" << temp.trans5() << "\t" << temp.mid5() << "\t"
				<< temp.cis3() << "\t" << temp.trans3() << "\t" << temp.mid3() << "\t"
				<< trueHelices[j]->freeEnergy(eStack, *alignedSeqs[i]) << "\t" << competingHelices[i]->at(k)->freeEnergy(eStack, *seqs[i]) << "\n";

			}

		}

	}
}

void Alignment::makeCompHelixCores(){
	vector<vector<int>* > bpCounts;
	//bpCounts[x][y] = number of times bp x-y appears in competing helices
	//across all sequences
	//could save space with triangular matrix, but not a huge concern at the moment

	//initialize bpCounts Matrix;
	for(unsigned int i = 0; i < alignedSeqs[0]->length(); i++){
		vector<int>* temp = new vector<int>();

		for(unsigned int j = 0; j < alignedSeqs[0]->length(); j++){
			temp->push_back(0);
		}
		bpCounts.push_back(temp);
	}

	//fill bpCounts Matrix:
	//loop through sequences: i
	for(unsigned int i = 0; i < competingHelices.size(); i++){
		//loop through competing helices: j
		for (unsigned int j = 0; j < competingHelices[i]->size(); j++){
			CompetingHelix temp = *(competingHelices[i]->at(j));
			int bp5 = temp.start5;
			int bp3 = temp.end3;
			while(bp5 <= temp.end5){
				int alignedBp5 = seq2AlignmentMap[i]->at(bp5);
				int alignedBp3 = seq2AlignmentMap[i]->at(bp3);

				bpCounts[alignedBp5]->at(alignedBp3) += 1;

				//next bp
				bp5++;
				bp3--;
			}
		}
	}


	//initialize core vectors:
	for(unsigned int i = 0; i <= alignedSeqs.size(); i++){
		cores.push_back(new vector<HelixCore*>());
	}

	//find cores
	for(int cutoff = (int)alignedSeqs.size(); cutoff >= ((int)alignedSeqs.size() * 3) / 4; cutoff--){
		for(int i = 0; i < (int)bpCounts.size(); i++){
			//starting from top right corner of bpCounts matrix
			for(int j = (int)bpCounts[i]->size() - 1; j > i; j--){

				if(bpCounts[i]->at(j) == cutoff){
					HelixCore * temp = new HelixCore(cutoff, i, j);
					//extend core:
					int extendI = i+1;
					int extendJ = j-1;
					while(bpCounts[extendI]->at(extendJ) == cutoff){
						temp->growIn();
						extendI++;
						extendJ--;
					}
					//figure out to what level counts should be reduced
					int outScore = 0; //bpCount outward of core
					int inScore = bpCounts[extendI]->at(extendJ); //bpCount inward of core
					if(i-1 >=0 && j+1 < (int)bpCounts[i]->size()){
						outScore = bpCounts[i-1]->at(j+1);
					}
					int newScore = max(inScore, outScore);

					//reduce bpCounts for core helix
					for(int k = 0; k < temp->getLength(); k++){
						bpCounts[i+k]->at(j-k) = newScore;
					}
					if(temp->getLength() > 0){
						cores[cutoff]->push_back(temp);
					}
					else{
						delete temp;
					}

				}

			}
		}
	}

	/*
	//print out core lengths for each cutoff:
	for(unsigned int i = 0; i < cores.size(); i++){
		if (cores[i]->size() == 0){
			//skip
		}
		else{
			cout << "Cutoff " << i << endl;
			double avg = 0;
			for(unsigned int j = 0; j < cores[i]->size(); j++){
				avg += cores[i]->at(j)->getLength();
				cout << cores[i]->at(j)->getLength() << endl;
			}

			avg = avg / cores[i]->size();
			cout << "avg = " << avg << endl;
		}
	}
	*/

	//free memory
	for(unsigned int i = 0; i < bpCounts.size(); i++){
		delete bpCounts[i];
	}
	bpCounts.clear();


}

//not finished yet...
bool Alignment::addHelixToGroups(vector<HelixGroup*> & existingGroups, CompetingHelix & newHelix, int newHelixSeqIndex)
{
	return false;
}


//not finished yet...
bool Alignment::isCompatible(HelixGroup & group, CompetingHelix & helix, int helixSeqIndex){
	//get helix core:
	int start5core = -1;
	int end5core = -1;
	int start3core = -1;
	int end3core = -1;

	for(unsigned int i = 0; i < group.getGroup().size(); i++){
		if (group.getGroup().at(i) != NULL){
			int start5current = seq2AlignmentMap[i]->at(group.getGroup().at(i)->start5);
			int end5current = seq2AlignmentMap[i]->at(group.getGroup().at(i)->end5);
			int start3current = seq2AlignmentMap[i]->at(group.getGroup().at(i)->start3);
			int end3current = seq2AlignmentMap[i]->at(group.getGroup().at(i)->end3);

			if(start5core == -1){
				start5core = start5current;
				end5core = end5current;
				start3core = start3current;
				end3core = end3current;
			}
			else{
				if(start5core < start5current){
					//shrink core
					start5core = start5current;
					end3core = end3current;
				}
				if(end5core > end5current){
					//shrink core
					end5core = end5current;

				}
			}
		}
	}

	return false;
}


void Alignment::shuffleAlignment()
{
	vector<int> unpairedColumns;
	//find unpaired columns:
	for(int i = 0; i < (int)alignedStruct.size(); i++){
		if(alignedStruct[i] == -1){
			unpairedColumns.push_back(i);
		}
	}
	//cout << "number of unpaired columns\n";

	srand((unsigned)time(0));
	//cout << "The value of RAND_MAX is " <<  RAND_MAX << endl;

	//shuffle list unpaired columns:
	for(unsigned int i = 0; i < unpairedColumns.size(); i++){
		unsigned int random_column = rand();

		random_column = i + random_column % (unpairedColumns.size()-i);

		//swap i with random_column
		int temp = unpairedColumns[i];
		unpairedColumns[i] = unpairedColumns[random_column];
		unpairedColumns[random_column] = temp;
	}

	//reorder columns:
	vector<string*> newAlignment;
	for(unsigned int i = 0; i < alignedSeqs.size(); i++){
		newAlignment.push_back(new string(*(alignedSeqs[i])));
	}

	int unpairedColsIndex = 0;
	for(unsigned int i = 0; i < alignedStruct.size(); i++){
		if(alignedStruct[i] == -1){
			for(unsigned int j = 0; j < alignedSeqs.size(); j++){
				newAlignment[j]->at(i) = alignedSeqs[j]->at(unpairedColumns[unpairedColsIndex]);
			}
			unpairedColsIndex++;
		}
	}

	//delete old alignedSeqs strings
	for (unsigned int i = 0; i < alignedSeqs.size(); i++){
		delete alignedSeqs[i];
	}
	alignedSeqs.clear();
	alignedSeqs = newAlignment;

	//clear old seqs information
	for(unsigned int i = 0; i < seqs.size(); i++){
		delete seqs[i];
	}
	seqs.clear();
	for(unsigned int i = 0; i < seq2AlignmentMap.size(); i++){
		delete seq2AlignmentMap[i];
	}
	seq2AlignmentMap.clear();
	for(unsigned int i = 0; i < alignment2SeqMap.size(); i++){
		delete alignment2SeqMap[i];
	}
	alignment2SeqMap.clear();

	//set seqs to reflect new shuffled alignment
	vector<int>* temp;
	vector<int>* temp2;
	string * seq;

	for (unsigned int i = 0; i < alignedSeqs.size(); i++){

		temp = new vector<int>();
		temp2 = new vector<int>();
		seq = new string();

		for (unsigned int j = 0; j < alignedSeqs[i]->length(); j++){
			if (alignedSeqs[i]->at(j) != '-'){
				temp->push_back(j);
				*seq = *seq + alignedSeqs[i]->at(j);
				temp2->push_back(temp->size() - 1);
			}
			else{
				temp2->push_back(-1);
			}
		}
		seqs.push_back(seq);
		seq2AlignmentMap.push_back(temp);
		alignment2SeqMap.push_back(temp2);

	}

	fillSeqStructs();
	labelHelices();

//	cout << alignedStructString() << endl;
//	cout << alignedlabelString() << endl;
//	cout << structString();
//	cout << labelString();
//
//	cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";

}

bool Alignment::helixIsClean(CompetingHelix & helix, int seqIndex){
	int pos5 = helix.start5;
	int pos3 = helix.end3;
	int length = helix.end5 - helix.start5 + 1;

	int alignedPos5, alignedPos3;
	for(int k = 0; k < length; k++){
		alignedPos5 = seq2AlignmentMap[seqIndex]->at(pos5+k);
		alignedPos3 = seq2AlignmentMap[seqIndex]->at(pos3-k);

		if(alignedStruct[alignedPos5] != -1 && alignedStruct[alignedPos3] != -1
				&& alignmentLabels[alignedPos5] == alignmentLabels[alignedPos3]){
			return false;
		}
	}
	return true;

}

void Alignment::nickTemp(Tree & tree, unsigned int randomizedTrials){
	double eStack[4][4][4][4];
	Utilities::ReadStack(eStack);



	for(unsigned int rand_index = 0; rand_index < randomizedTrials; rand_index++){
		ShuffledAlignment shuffled(*this, 1);

		shuffled.alignmentGStats(eStack);
		for(unsigned int seq_index = 0; seq_index < shuffled.competingHelices.size(); seq_index++){
			for(unsigned int helix_index = 0; helix_index < shuffled.competingHelices[seq_index]->size(); helix_index++){
				if(!shuffled.StatsMatrix[seq_index]->at(0)->at(helix_index)->isZero() && shuffled.helixIsClean(*shuffled.competingHelices[seq_index]->at(helix_index), seq_index))
					cout << shuffled.logLikelihood(tree, *shuffled.competingHelices[seq_index]->at(helix_index), seq_index) << endl;
			}

		}

	}

}

void Alignment::calculatePValues(Tree & tree, string & treeFile, unsigned int randomizedTrials){
	if(printHeaders){
		cout << "p-value\tLog Likelihood\tsequence index\ttrue helix\tcompeting helix\t"
		<< "cis3\tcis5\ttrans3\ttrans5\tmid3\tmid5\t"
		<< "midpoint\tstructure\t"
		<< "%consensusBP\t%noGapConsensusBP\t%gaps\tseqCons\tCovariance"
		<< endl;
	}

	assert(randomizedTrials > 0);//sanity check

	double eStack[4][4][4][4];
	Utilities::ReadStack(eStack);


	//loop through true helices, generating random alignments
	for(int i = 0; i < helixNumber; i++){

		//cout << "True Helix " << i << endl;
		//cout << trueHelices[i]->printHelix(alignedStruct.size()) << endl;

		//figure out which competing Helices are relevant (i.e. which ones compete with
		//true helix i
		vector<pair<CompetingHelix*, int> > relevantHelices;
		vector<double> logLikes;
		vector<double> pvalues;
		vector<StatsWrapper*> statsWrappers;

		for(unsigned int j = 0; j < StatsMatrix.size(); j++){
			for(unsigned int k = 0; k < StatsMatrix[j]->at(i)->size(); k++){
				if(!StatsMatrix[j]->at(i)->at(k)->isZero() && helixIsClean(*competingHelices[j]->at(k), j)){
					relevantHelices.push_back(pair<CompetingHelix*, int>(competingHelices[j]->at(k), j));
					statsWrappers.push_back(StatsMatrix[j]->at(i)->at(k));
					logLikes.push_back(logLikelihood(tree, *competingHelices[j]->at(k), j));
					pvalues.push_back(0);
				}
			}
		}
		//cout << "realigning...\n";
		ShuffledAlignment * realigned = new ShuffledAlignment(*this, i, treeFile);
		//perform randomized trials:
		vector<double> nullDistLogLikes;
		for(unsigned int rand_index = 0; rand_index < randomizedTrials; rand_index++){
			//cerr << "helix " << i << " trial: " << rand_index << endl;

			if(rand_index % 100 == 0 && rand_index != 0){
				delete realigned;
				//cout << "realigning...\n";
				realigned = new ShuffledAlignment(*this, i, treeFile);
			}
			ShuffledAlignment shuffled(*realigned, 0); //zero because the realigned alignment has only the one true helix



			shuffled.alignmentGStats(eStack);
			for(unsigned int seq_index = 0; seq_index < shuffled.competingHelices.size(); seq_index++){
				for(unsigned int helix_index = 0; helix_index < shuffled.competingHelices[seq_index]->size(); helix_index++){
					if(!shuffled.StatsMatrix[seq_index]->at(0)->at(helix_index)->isZero() && shuffled.helixIsClean(*shuffled.competingHelices[seq_index]->at(helix_index), seq_index))
						nullDistLogLikes.push_back(shuffled.logLikelihood(tree, *shuffled.competingHelices[seq_index]->at(helix_index), seq_index));
				}

			}

			sort(nullDistLogLikes.begin(), nullDistLogLikes.end());

			for(unsigned int k = 0; k < pvalues.size(); k++){
				if(nullDistLogLikes.size() > 0){
					pvalues[k] += 1 - Utilities::lowerBound(nullDistLogLikes, 0, nullDistLogLikes.size(), logLikes[k])/(double) nullDistLogLikes.size();
				}
				else{
					//implicit: pvalues[k] += 0;
				}
			}

			nullDistLogLikes.clear();
		}

		delete realigned;

		for(unsigned int k = 0; k < pvalues.size(); k++){
			pvalues[k] = pvalues[k] / randomizedTrials;
			cout << pvalues[k] << "\t" << logLikes[k] << "\t" << relevantHelices[k].second <<"\t"<< i << "\t" << k << "\t";
			cout << statsWrappers[k]->cis3() << "\t" << statsWrappers[k]->cis5() << "\t";
			cout << statsWrappers[k]->trans3() << "\t" << statsWrappers[k]->trans5() << "\t";
			cout << statsWrappers[k]->mid3() << "\t" << statsWrappers[k]->mid5() << "\t";
			cout << relevantHelices[k].first->midpoint() << "\t";

			cout << relevantHelices[k].first->printAlignedHelix(*seq2AlignmentMap[relevantHelices[k].second], alignedStruct.size()) << "\t";
			cout << consensusBPPercent(*relevantHelices[k].first, relevantHelices[k].second) << "\t";
			cout << noGapCBPPercent(*relevantHelices[k].first, relevantHelices[k].second) << "\t";
			cout << getGapFraction(*relevantHelices[k].first, relevantHelices[k].second) << "\t";
			cout << getSeqCons(*relevantHelices[k].first, relevantHelices[k].second) << "\t";
			cout << getCovariance(*relevantHelices[k].first, relevantHelices[k].second);
			cout << endl;


		}

	}
}

void Alignment::calculatePValues(Tree & tree, unsigned int randomizedTrials){

	if(printHeaders){
		cout << "p-value\tLog Likelihood\tsequence index\ttrue helix\tcompeting helix\t"
		<< "cis3\tcis5\ttrans3\ttrans5\tmid3\tmid5\t"
		<< "midpoint\tstructure\n";
	}

	assert(randomizedTrials > 0);//sanity check

	double eStack[4][4][4][4];
	Utilities::ReadStack(eStack);


	//loop through true helices, generating random alignments
	for(int i = 0; i < helixNumber; i++){

		//cout << "True Helix " << i << endl;
		//cout << trueHelices[i]->printHelix(alignedStruct.size()) << endl;

		//figure out which competing Helices are relevant (i.e. which ones compete with
		//true helix i
		vector<pair<CompetingHelix*, int> > relevantHelices;
		vector<double> logLikes;
		vector<double> pvalues;
		vector<StatsWrapper*> statsWrappers;

		for(unsigned int j = 0; j < StatsMatrix.size(); j++){
			for(unsigned int k = 0; k < StatsMatrix[j]->at(i)->size(); k++){
				if(!StatsMatrix[j]->at(i)->at(k)->isZero() && helixIsClean(*competingHelices[j]->at(k), j)){
					relevantHelices.push_back(pair<CompetingHelix*, int>(competingHelices[j]->at(k), j));
					statsWrappers.push_back(StatsMatrix[j]->at(i)->at(k));
					logLikes.push_back(logLikelihood(tree, *competingHelices[j]->at(k), j));
					pvalues.push_back(0);
				}
			}
		}

		//perform randomized trials:
		vector<double> nullDistLogLikes;
		for(unsigned int rand_index = 0; rand_index < randomizedTrials; rand_index++){
			//cerr << "helix " << i << " trial: " << rand_index << endl;
			ShuffledAlignment shuffled(*this, i);

			shuffled.alignmentGStats(eStack);
			for(unsigned int seq_index = 0; seq_index < shuffled.competingHelices.size(); seq_index++){
				for(unsigned int helix_index = 0; helix_index < shuffled.competingHelices[seq_index]->size(); helix_index++){
					if(!shuffled.StatsMatrix[seq_index]->at(0)->at(helix_index)->isZero() && shuffled.helixIsClean(*shuffled.competingHelices[seq_index]->at(helix_index), seq_index))
						nullDistLogLikes.push_back(shuffled.logLikelihood(tree, *shuffled.competingHelices[seq_index]->at(helix_index), seq_index));
				}

			}

			sort(nullDistLogLikes.begin(), nullDistLogLikes.end());

			for(unsigned int k = 0; k < pvalues.size(); k++){
				pvalues[k] += 1 - Utilities::lowerBound(nullDistLogLikes, 0, nullDistLogLikes.size(), logLikes[k])/(double) nullDistLogLikes.size();
			}

			nullDistLogLikes.clear();
		}

		for(unsigned int k = 0; k < pvalues.size(); k++){
			pvalues[k] = pvalues[k] / randomizedTrials;
			cout << pvalues[k] << "\t" << logLikes[k] << "\t" << relevantHelices[k].second <<"\t"<< i << "\t" << k << "\t";
			cout << statsWrappers[k]->cis3() << "\t" << statsWrappers[k]->cis5() << "\t";
			cout << statsWrappers[k]->trans3() << "\t" << statsWrappers[k]->trans5() << "\t";
			cout << statsWrappers[k]->mid3() << "\t" << statsWrappers[k]->mid5() << "\t";
			cout << relevantHelices[k].first->midpoint() << "\t";
			if(pvalues[k] < 0.05){
				cout << relevantHelices[k].first->printAlignedHelix(*seq2AlignmentMap[relevantHelices[k].second], alignedStruct.size());
			}
			else{
				cout << "-";
			}
			cout << endl;


		}

		//benjamini-hochberg correction:
//		vector<double> adjustedPvalues(pvalues);
//
//		for(unsigned int k = 0; k < adjustedPvalues.size();k++){
//			int rank = adjustedPvalues.size();
//			for(unsigned int l = 0; l < adjustedPvalues.size();l++){
//				if (pvalues[k] < pvalues[l]){
//					rank--;
//				}
//			}
//			adjustedPvalues[k] *= adjustedPvalues.size() / (double)rank;
//
//			//cout << rank << "\t" << logLikes[k] << "\t" << pvalues[k] << "\t" << adjustedPvalues[k] << "\t" << meanValidBPs(*relevantHelices[k].first, relevantHelices[k].second) << endl;
//			if(adjustedPvalues[k] < 0.05){
//				cout << "Log Likelihood ratio: " << logLikes[k] << endl;
//				cout << "P-value: " << pvalues[k] << endl;
//				cout << "Rank: " << rank << endl;
//				cout << "Adjusted P-value: " << adjustedPvalues[k] << endl;
//				cout << "Mean valid BPs per column: " << meanValidBPs(*relevantHelices[k].first, relevantHelices[k].second) << endl;
//
//				cout << relevantHelices[k].first->printAlignedHelix(*seq2AlignmentMap[relevantHelices[k].second], alignedStruct.size()) << endl;
//				cout << printAlignmentAtCompetingHelix(*relevantHelices[k].first, relevantHelices[k].second);
//
//				cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
//			}
//
//		}

	}
}
//
//void Alignment::calculateCorePValues(Tree & tree, unsigned int randomizedTrials){
//
//	assert(randomizedTrials > 0);//sanity check
//
//	srand((unsigned)time(0));
//
//	double eStack[4][4][4][4];
//	Utilities::ReadStack(eStack);
//
//
//	//loop through true helices, generating random alignments
//	for(unsigned int i = 0; i < trueHelices.size(); i++){
//
//		cout << "True Helix " << i << endl;
//		cout << trueHelices[i]->printHelix(alignedStruct.size()) << endl;
//
//		//figure out which competing Helices are relevant (i.e. which ones compete with
//		//true helix i
//		vector<pair<CompetingHelix*, int> > relevantHelices;
//		vector<double> logLikes;
//		vector<double> pvalues;
//
//		for(unsigned int j = 0; j < StatsMatrix.size(); j++){
//			for(unsigned int k = 0; k < StatsMatrix[j]->at(i)->size(); k++){
//				if(!StatsMatrix[j]->at(i)->at(k)->isZero()){
//					relevantHelices.push_back(pair<CompetingHelix*, int>(competingHelices[j]->at(k), j));
//					logLikes.push_back(logLikelihood(tree, *competingHelices[j]->at(k), j));
//					pvalues.push_back(0);
//				}
//			}
//		}
//
//		vector<double> nullDistLogLikes;
//		for(unsigned int rand_index = 0; rand_index < randomizedTrials; rand_index++){
//			ShuffledAlignment shuffled(*this, i);
//
//			shuffled.alignmentGStats(eStack);
//			for(unsigned int seq_index = 0; seq_index < shuffled.competingHelices.size(); seq_index++){
//				for(unsigned int helix_index = 0; helix_index < shuffled.competingHelices[seq_index]->size(); helix_index++){
//					if(!shuffled.StatsMatrix[seq_index]->at(0)->at(helix_index)->isZero())
//						nullDistLogLikes.push_back(shuffled.logLikelihood(tree, *shuffled.competingHelices[seq_index]->at(helix_index), seq_index));
//				}
//
//			}
//
//			sort(nullDistLogLikes.begin(), nullDistLogLikes.end());
//
//			for(unsigned int k = 0; k < pvalues.size(); k++){
//				pvalues[k] += 1 - Utilities::lowerBound(nullDistLogLikes, 0, nullDistLogLikes.size(), logLikes[k])/(double) nullDistLogLikes.size();
//			}
//
//			nullDistLogLikes.clear();
//		}
//
//		for(unsigned int k = 0; k < pvalues.size(); k++){
//			pvalues[k] = pvalues[k] / randomizedTrials;
////			cout << pvalues[k] << endl;
//		}
//
//		vector<double> adjustedPvalues(pvalues);
//
//		for(unsigned int k = 0; k < adjustedPvalues.size();k++){
//			int rank = adjustedPvalues.size();
//			for(unsigned int l = 0; l < adjustedPvalues.size();l++){
//				if (pvalues[k] < pvalues[l]){
//					rank--;
//				}
//			}
//			adjustedPvalues[k] *= adjustedPvalues.size() / (double)rank;
//
//			cout << rank << "\t" << logLikes[k] << "\t" << pvalues[k] << "\t" << adjustedPvalues[k] << "\t" << meanValidBPs(*relevantHelices[k].first, relevantHelices[k].second) << endl;
//			if(adjustedPvalues[k] < 0.05){
//				cout << "Log Likelihood ratio: " << logLikes[k] << endl;
//				cout << "P-value: " << pvalues[k] << endl;
//				cout << "Rank: " << rank << endl;
//				cout << "Adjusted P-value: " << adjustedPvalues[k] << endl;
//
//				cout << relevantHelices[k].first->printAlignedHelix(*seq2AlignmentMap[relevantHelices[k].second], alignedStruct.size()) << endl;
//
//
//				cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
//			}
//
//		}
//
//	}
//}

double Alignment::meanValidBPs(CompetingHelix & h, int seqIndex){
	int pos5 = h.start5;
	int pos3 = h.end3;
	int length = h.end5 - pos5 + 1;

	int bp_count = 0;
	for (int k = 0; k < length; k++){
		int alignedPos5 = seq2AlignmentMap[seqIndex]->at(pos5+k);
		int alignedPos3 = seq2AlignmentMap[seqIndex]->at(pos3-k);
		for(unsigned int i = 0; i < alignedSeqs.size(); i++){
			char a = alignedSeqs[i]->at(alignedPos5);
			char b = alignedSeqs[i]->at(alignedPos3);
			if (Utilities::validBP(a,b)){
				bp_count++;
			}
		}
	}
	return bp_count / (double)(alignedSeqs.size() * length);
}

double Alignment::logLikelihood(Tree & treeRoot, CompetingHelix & h, int seqIndex){

	double logLike = 0;
	int pos5 = h.start5;
	int pos3 = h.end3;

	int length = h.end5 - pos5 + 1;
	//convert to alignment positions:
	for (int k = 0; k < length; k++){
		int alignedPos5 = seq2AlignmentMap[seqIndex]->at(pos5+k);
		int alignedPos3 = seq2AlignmentMap[seqIndex]->at(pos3-k);
		logLike += logLikelihood(alignedPos5 , alignedPos3, treeRoot);

	}

	return (logLike)/length;
}

string Alignment::printAlignmentAtCompetingHelix(CompetingHelix & h, int seqIndex){
	int alignedPos5 = seq2AlignmentMap[seqIndex]->at(h.start5);
	int alignedPos3 = seq2AlignmentMap[seqIndex]->at(h.end3);

	int alignedLength = seq2AlignmentMap[seqIndex]->at(h.end5) - alignedPos5 +1;

	string output = "";

	for(unsigned int j = 0; j < alignedSeqs.size(); j++){
		for(int i = 0; i < alignedLength; i++){
			output += alignedSeqs[j]->at(alignedPos5+i);
		}
		output+= "...";
		for(int i = alignedLength-1; i >=0; i--){

			output += alignedSeqs[j]->at(alignedPos3-i);
		}
		output+= "\n";
	}
	for(int i = 0; i < alignedLength; i++){
		if(alignment2SeqMap[seqIndex]->at(alignedPos5+i) != -1){
			output+= "(";
		}
		else{
			output+= "-";
		}

	}
	output+= "...";
	for(int i = alignedLength-1; i >=0; i--){
		if(alignment2SeqMap[seqIndex]->at(alignedPos3-i) != -1){
			output+= ")";
		}
		else{
			output+= "-";
		}
	}
	output+= "\n";
	return output;
}

/**
 * Current implementation is ultra-super inefficient - I needed to code it quickly...
 */
double Alignment::calculatePercentOfCompetingHelicesInCores(){

	int competingHelixCount = 0;
	int inCores = 0;

	bool found = false;

	for(unsigned int i = 0; i < competingHelices.size(); i++){

		for(unsigned int j = 0; j < competingHelices[i]->size(); j++){
			competingHelixCount++;
			found = false;

			int pos5 = competingHelices[i]->at(j)->start5;
			int pos3 = competingHelices[i]->at(j)->end3;

			int length = competingHelices[i]->at(j)->end5 - pos5 + 1;
			for(int k = 0; k < length && !found; k++){
				int alignedPos5 = seq2AlignmentMap[i]->at(pos5+k);
				int alignedPos3 = seq2AlignmentMap[i]->at(pos3-k);

				for(unsigned int l = 0; l < cores.size() && !found; l++){
					for(unsigned int m = 0; m < cores[l]->size() && !found; m++){
						int core5 = cores[l]->at(m)->getPos5();
						int core3 = cores[l]->at(m)->getPos3();
						int coreLength = cores[l]->at(m)->getLength();

						for(int n = 0; n < coreLength && !found; n++){
							if(core5+n == alignedPos5 && core3-n == alignedPos3){
								found = true;
								inCores++;
							}
						}

					}
				}
			}

		}
	}
	cout << inCores << endl;
	cout << competingHelixCount << endl;

	return inCores/ (double)competingHelixCount;

}

double Alignment::logLikelihood(int pos5, int pos3, Tree & tree){

	double paired = felsDoubles[pos5]->at(pos3);
	if(paired > 0){
		paired = log(tree.calcFelsDouble(*this, pos5, pos3))/log(2);
		assert(paired <= 0);
		felsDoubles[pos5]->at(pos3) = paired;
	}

	double pos5Single = felsSingles[pos5];
	if(pos5Single > 0){
		pos5Single = log(tree.calcFelsSingle(*this, pos5))/log(2);
		assert(pos5Single <= 0);
		felsSingles[pos5] = pos5Single;
	}

	double pos3Single = felsSingles[pos3];
	if(pos3Single > 0){
		pos3Single = log(tree.calcFelsSingle(*this, pos3)) / log(2);
		assert(pos3Single <= 0);
		felsSingles[pos3] = pos3Single;
	}

	return paired - pos5Single - pos3Single;
}

double Alignment::logPairedLikelihood(int pos5, int pos3, Tree & tree){
	double paired = felsDoubles[pos5]->at(pos3);
	if(paired > 0){
		paired = log(tree.calcFelsDouble(*this, pos5, pos3))/log(2);
		assert(paired <= 0);
		felsDoubles[pos5]->at(pos3) = paired;
	}
	return paired;
}

double Alignment::logUnpairedLikelihood(int pos, Tree & tree){
	double posSingle = felsSingles[pos];
	if(posSingle > 0){
		posSingle = log(tree.calcFelsSingle(*this, pos))/log(2);
		assert(posSingle <= 0);
		felsSingles[pos] = posSingle;
	}

	return posSingle;
}

vector<string> Alignment::realignTcoffee(Tree & tree) const{
	//write fasta file:
	string fastaFilenameTemplate = "/tmp/fasta.XXXXXX";
	char fastaFilename[fastaFilenameTemplate.size()+1];

	strcpy(fastaFilename, fastaFilenameTemplate.c_str());

	int unique_fd;
	FILE * unique_file;

	if ((unique_fd = mkstemp(fastaFilename)) == -1 ){
		cerr << "Could not create temporary clustal file for shuffling\n";
		exit(-1);
	}

	if((unique_file = fdopen(unique_fd, "w")) == NULL){
		cerr << "Could not open temporary clustal file for shuffling\n";
		exit(-1);
	}

	fputs(FastaFormat(false).c_str(), unique_file);
	fclose(unique_file);

	//make a temp file for output:
	string outTemplate = "/tmp/realign_out.XXXXXX";
	char outFilename[outTemplate.size()+1];
	strcpy(outFilename, outTemplate.c_str());

	if ((unique_fd = mkstemp(outFilename)) == -1 ){
		cerr << "Could not create temporary out file for realignment\n";
		exit(-1);
	}

	if((unique_file = fdopen(unique_fd, "w")) == NULL){
		cerr << "Could not open temporary out file for realignment\n";
		exit(-1);
	}
	string temp_contents = "TEMP!\n"; // this will be replaced by TCOFFEE output
	fputs(temp_contents.c_str(), unique_file);
	fclose(unique_file);

	//make a temp file for output:
	string treeTemplate = "/tmp/tree_out.XXXXXX";
	char treeFilename[treeTemplate.size()+1];
	strcpy(treeFilename, treeTemplate.c_str());

	if ((unique_fd = mkstemp(treeFilename)) == -1 ){
		cerr << "Could not create temporary tree file for realignment\n";
		exit(-1);
	}

	if((unique_file = fdopen(unique_fd, "w")) == NULL){
		cerr << "Could not open temporary out file for realignment\n";
		exit(-1);
	}
	//write tree to file
	fputs(tree.newickString().c_str(), unique_file);
	fclose(unique_file);

	//run realignment:
	string command = TCOFFEE_LOC + " " + string(fastaFilename) + " -outfile " + string(outFilename) + " -usetree " + string(treeFilename) + " -output fasta -n_core 1 -outorder input -quiet";
	if(system(command.c_str()) != 0){
		cerr << "Error: t_coffee encountered an error\n";
		exit(-1);
	}

	ifstream realignFile(outFilename);

	if(!realignFile.is_open()){
		cerr << "Error: cannot open realignment fasta file " << string(outFilename) << endl;
		exit(-1);
	}

	string line;
	vector<string> realignment;
	vector<string> r_names;
	string seq;
	bool found_seq = false;

	while(realignFile.good()){
		getline(realignFile, line);
		Utilities::TrimSpaces(line);
		//skip blank lines
		if(line.length() == 0){
			continue;
		}
		if(line.at(0) == '>'){
			if(found_seq){
				//cout << "Sequence Found: " << *seq << endl;
				realignment.push_back(seq);
			}
			line = line.substr(1);
			r_names.push_back(line);
			found_seq = true;
			seq = "";
		}
		else{
			seq = seq + line;
		}
	}

	realignFile.close();

	realignment.push_back(seq);


	assert(realignment.size() == r_names.size());
	for(unsigned int j = 0; j < r_names.size(); j++){

		//		cout << realignment[j] << endl;

		//		cout  << r_names[j] << "\t" << *seqNames[j] << endl;
		assert(r_names[j].compare(*seqNames[j]) == 0);
		assert(realignment[j].length() == realignment[0].length());
	}



	//clean up files
	command = "rm -f " + string(outFilename) + " " + string(fastaFilename) + " " + string(treeFilename) ;
	if(system(command.c_str()) != 0){
		cerr << "Error cleaning up temporary shuffle output file\n";
		exit(-1);
	}

	return realignment;
}

vector<string> Alignment::realignInterval(int begin, int end, string & treeFilename) const{

	assert(begin < end);

	//populate interval vector with empty lists

//	vector<list<char> > interval;
//	for(unsigned int j = 0; j < alignedSeqs.size(); j++){
//		interval.push_back(list<char>);
//	}

	vector<string> startingAlignment;
	for(unsigned int j = 0; j < alignedSeqs.size(); j++){
		startingAlignment.push_back(alignedSeqs[j]->substr(begin, end-begin));
	}

	string tempFilename = alignmentName + "." + TEMP_FASTA_FILENAME;
	string tempOutputFilename = alignmentName + "." + TEMP_OUTPUT_FASTA_FILENAME;

	tempFilename = Utilities::checkForConflicts(tempFilename);
	tempOutputFilename = Utilities::checkForConflicts(tempOutputFilename);

	ofstream tempFasta;
	tempFasta.open(tempFilename.c_str());
	assert(tempFasta.good());

	for(unsigned int j = 0; j < startingAlignment.size(); j++){
		tempFasta << ">" << *seqNames[j] << endl;
		tempFasta << startingAlignment[j] << endl;
	}

	tempFasta.close();

	string command;
	//command = "java -version";
	//system(command.c_str());

	command = JAVA_LOC + " -cp " + REALIGNER_LOCATION + " Realigner " + tempFilename + " "
	+ treeFilename + " > " + tempOutputFilename;

	int out_val = system(command.c_str());

	//cerr << "Java Out: " << out_val << endl;
	assert(out_val == 0);

	ifstream realignFile(tempOutputFilename.c_str());

	if(!realignFile.is_open()){
		cerr << "Error: cannot open realignment fasta file " << tempOutputFilename << endl;
		exit(-1);
	}

	string line;
	vector<string> realignment;
	vector<string> r_names;
	string seq;
	bool found_seq = false;

	while(realignFile.good()){
		getline(realignFile, line);
		Utilities::TrimSpaces(line);
		//skip blank lines
		if(line.length() == 0){
			continue;
		}
		if(line.at(0) == '>'){
			if(found_seq){
				//cout << "Sequence Found: " << *seq << endl;
				realignment.push_back(seq);
			}
			line = line.substr(1);
			r_names.push_back(line);
			found_seq = true;
			seq = "";
		}
		else{
			seq = seq + line;
		}
	}

	realignFile.close();

	realignment.push_back(seq);


	assert(realignment.size() == r_names.size());
	for(unsigned int j = 0; j < r_names.size(); j++){

//		cout << realignment[j] << endl;

//		cout  << r_names[j] << "\t" << *seqNames[j] << endl;
		assert(r_names[j].compare(*seqNames[j]) == 0);
		assert(realignment[j].length() == realignment[0].length());
	}

	//clean up temp files:
	command = "rm -f " + tempFilename;
	system(command.c_str());
	command = "rm -f " + tempOutputFilename;
	system(command.c_str());

	return realignment;

}

void Alignment::printTrueHelixStats(Tree & tree){

	if(Alignment::printHeaders){
		cout << "%ConsensusBP\t%ConsensusBP(gap pairs omitted)\tgap fraction\tLogLikelihood\tLength\tSeqConservation\tCovariance\tStructure\n";
	}

	for(int i = 0; i < helixNumber; i++){

		double consensusBP = 0;
		double totalBP = 0;
		double noGapCBP = 0;
		double noGapTBP = 0;
		double gapCount = 0;
		double seqCons = 0;
		double covariance = 0;

		double loglike =  0;
		int length = 0;

		for(unsigned int j = 0; j < alignmentLabels.size(); j++){

			//found a position in helix i
			if(alignmentLabels[j] == i && alignedStruct[j] > (int)j){
				//find pairing partner:
				int pp = alignedStruct[j];
				assert(pp >= 0);

				loglike += 	logLikelihood(j, pp, tree);
				length++;

				for(unsigned int k = 0; k < alignedSeqs.size(); k++){
					char a, b;
					a = alignedSeqs[k]->at(j);
					b = alignedSeqs[k]->at(pp);

					if(Utilities::validBP(a,b)){
						consensusBP++;
					}
					totalBP++;

					if(a != '-' && b != '-'){
						if(Utilities::validBP(a,b)){
							noGapCBP++;
						}
						noGapTBP++;
					}
					else{
						if( a == '-'){
							gapCount++;
						}
						if(b == '-'){
							gapCount++;
						}
					}

					seqCons += getSeqCons(j);
					seqCons += getSeqCons(pp);
					covariance += getCovariance(j,pp);
				}

			}
		}

		seqCons = seqCons / (double)(totalBP * 2);
		covariance = covariance / totalBP;

		cout << consensusBP / (double) totalBP << "\t" << noGapCBP / (double) noGapTBP << "\t" << gapCount / (double)(totalBP*2) << "\t"
		<< loglike / length << "\t" << length << "\t" << seqCons << "\t" << covariance <<"\t" << TrueHelixString(i)
		<< endl;
	}
}

double Alignment::consensusBPPercent(int pos5, int pos3){

	if(consensusBP == NULL){
		consensusBP = new UTMatrix(alignedStruct.size());

		for(unsigned int i = 0; i < alignedStruct.size(); i++){
			for(unsigned int j = i; j < alignedStruct.size(); j++){
				consensusBP->set(i, j, -1);
			}
		}
	}

	if(consensusBP->at(pos5, pos3) < 0){
		int total = 0;
		int cbp = 0;
		double cbp_percent;

		for(unsigned int i = 0; i < alignedSeqs.size(); i++){
			char a,b;
			a = alignedSeqs[i]->at(pos5);
			b = alignedSeqs[i]->at(pos3);

			total++;
			if(Utilities::validBP(a, b)){
				cbp++;
			}
		}
		assert(total > 0);
		cbp_percent = cbp / (double) total;
		consensusBP->set(pos5, pos3, cbp_percent);
		return cbp_percent;
	}
	else{
		return consensusBP->at(pos5, pos3);
	}

}
//double Alignment::noGapCBPPercent(int pos5, int pos3){
//
//	if(noGapCBP == NULL){
//		noGapCBP  = new UTMatrix(alignedStruct.size());
//
//		for(unsigned int i = 0; i < alignedStruct.size(); i++){
//			for(unsigned int j = i; j < alignedStruct.size(); j++){
//				noGapCBP->set(i, j, -1);
//			}
//		}
//	}
//
//	if(noGapCBP->at(pos5, pos3) < 0){
//		int total = 0;
//		int cbp = 0;
//		double cbp_percent;
//
//		for(unsigned int i = 0; i < alignedSeqs.size(); i++){
//			char a,b;
//			a = alignedSeqs[i]->at(pos5);
//			b = alignedSeqs[i]->at(pos3);
//
//			if(a != '-' && b != '-'){
//				total++;
//				if(validBP(a, b)){
//					cbp++;
//				}
//			}
//		}
//
//		assert(total > 0);
//		cbp_percent = cbp / (double) total;
//		noGapCBP->set(pos5, pos3, cbp_percent);
//		return cbp_percent;
//	}
//	else{
//		return noGapCBP->at(pos5, pos3);
//	}
//
//}

double Alignment::consensusBPPercent(CompetingHelix & helix, int seqIndex){
	int length = helix.end5 - helix.start5 + 1;
	int pos5;
	int pos3;

	double cbp = 0;

	for(int i = 0; i < length; i++){
		pos5 = seq2AlignmentMap[seqIndex]->at(helix.start5 + i);
		pos3 = seq2AlignmentMap[seqIndex]->at(helix.end3 - i);

		cbp += consensusBPPercent(pos5, pos3);
	}

	return cbp / length;
}

double Alignment::noGapCBPPercent(CompetingHelix & helix, int seqIndex){
	int length = helix.end5 - helix.start5 + 1;
	int pos5;
	int pos3;

	int total = 0;
	int cbp = 0;


	for(int i = 0; i < length; i++){
		pos5 = seq2AlignmentMap[seqIndex]->at(helix.start5 + i);
		pos3 = seq2AlignmentMap[seqIndex]->at(helix.end3 - i);

		for(unsigned int j = 0; j < alignedSeqs.size(); j++){
			char a,b;
			a = alignedSeqs[j]->at(pos5);
			b = alignedSeqs[j]->at(pos3);

			if(a != '-' && b != '-'){
				total++;
				if(Utilities::validBP(a, b)){
					cbp++;
				}
			}
		}

	}

	assert(total > 0);
	return cbp / (double) total;
}

double Alignment::getGapFraction(unsigned int pos){
	if (gapFraction == NULL){
		gapFraction = new vector<double>();
		for(unsigned int i = 0; i < alignedStruct.size(); i++){
			gapFraction->push_back(-1);
		}
	}
	if(gapFraction->at(pos) < 0){
		int gaps = 0;
		int total = 0;
		double gf;
		for(unsigned int i = 0; i < alignedSeqs.size(); i++){
			char a = alignedSeqs[i]->at(pos);
			if(a == '-'){
				gaps++;
			}
			total++;
		}
		gf = gaps / (double)total;
		(*gapFraction)[pos] = gf;
		return gf;
	}
	else{
		return gapFraction->at(pos);
	}
}

double Alignment::getGapFraction(CompetingHelix & helix, int seqIndex){
	int length = helix.end5 - helix.start5 + 1;
	int pos5;
	int pos3;

	double gf = 0;

	for(int i = 0; i < length; i++){
		pos5 = seq2AlignmentMap[seqIndex]->at(helix.start5 + i);
		pos3 = seq2AlignmentMap[seqIndex]->at(helix.end3 - i);

		gf += getGapFraction(pos5);
		gf += getGapFraction(pos3);
	}

	return gf / (length*2);
}

double Alignment::getSeqCons(unsigned int pos){
	if (seqCons == NULL){
		seqCons = new vector<double>();
		for(unsigned int i = 0; i < alignedStruct.size(); i++){
			seqCons->push_back(-1);
		}
	}
	if(seqCons->at(pos) < 0){

		int total = 0;
		int matching = 0;
		for(unsigned int i = 0; i < alignedSeqs.size() - 1; i++){
			for(unsigned int j = i + 1; j < alignedSeqs.size(); j++){
				if(toupper(alignedSeqs[i]->at(pos)) == toupper(alignedSeqs[j]->at(pos))){
					matching++;
				}
				total++;
			}
		}
		(*seqCons)[pos] = matching / (double) total;
		//(*seqCons)[pos] = seqConservationAtPosition(pos);
	}
	return seqCons->at(pos);

}

double Alignment::getSeqCons(CompetingHelix & helix, int seqIndex){
	int length = helix.end5 - helix.start5 + 1;
	int pos5;
	int pos3;

	double sc = 0;

	for(int i = 0; i < length; i++){
		pos5 = seq2AlignmentMap[seqIndex]->at(helix.start5 + i);
		pos3 = seq2AlignmentMap[seqIndex]->at(helix.end3 - i);

		sc += getSeqCons(pos5);
		sc += getSeqCons(pos3);
	}

	return sc / (length*2);
}

string Alignment::TrueHelixString(int th_index){
	string out = "";
	for(unsigned int i = 0; i < alignmentLabels.size(); i++){
		if(alignmentLabels[i] != th_index){
			out += ".";
		}
		else if(alignedStruct[i] > (int)i){
			out += "(";
		}
		else{
			out += ")";
		}
	}
	return out;
}

double Alignment::getCovariance(int pos5, int pos3){
	int denom = 0;
	int covariance = 0;
	for(unsigned int i = 0; i < alignedSeqs.size() - 1; i++){

		for(unsigned int j = i+1; j < alignedSeqs.size(); j++){
			int hamming = 0;
			char ia, ib, ja, jb;
			ia = toupper(alignedSeqs[i]->at(pos5));
			ib = toupper(alignedSeqs[i]->at(pos3));

			ja = toupper(alignedSeqs[j]->at(pos5));
			jb = toupper(alignedSeqs[j]->at(pos3));
			if(ia == ja){
				hamming++;
			}
			if(ib == jb){
				hamming++;
			}

			if(hamming > 0){
				if(Utilities::validBP(ia, ib) && Utilities::validBP(ja,jb)){
					covariance += hamming;
				}
				else{
					covariance -= hamming;
				}
			}
			denom++;
		}
	}
	return covariance / (double)denom;
}

double Alignment::getCovariance(CompetingHelix & helix, int seqIndex){
	int length = helix.end5 - helix.start5 + 1;
		int pos5;
		int pos3;

		double cov = 0;

		for(int i = 0; i < length; i++){
			pos5 = seq2AlignmentMap[seqIndex]->at(helix.start5 + i);
			pos3 = seq2AlignmentMap[seqIndex]->at(helix.end3 - i);

			cov += getCovariance(pos5, pos3);
		}

		return cov / length;
}



string Alignment::ClustalWFormat(){
	stringstream out;
	out << "CLUSTAL W\n\n";

	for(unsigned int i = 0; i < alignedSeqs.size();i++){
		stringstream temp;
		temp << i;
		assert(temp.str().length() <= 10);
		for(unsigned int j = temp.str().length(); j < 11; j++){
			temp << " ";
		}
		temp << *alignedSeqs[i];
		out << temp.str() << endl;
	}
	//print line with *s
	for(unsigned int j = 0; j < 11; j++){
		out << " ";
	}
	for(unsigned int j = 0; j < alignedStruct.size(); j++){
		bool conserved = true;
		char first_char = alignedSeqs[0]->at(j);
		for(unsigned int i = 1; i < alignedSeqs.size();i++){
			if(alignedSeqs[i]->at(j) != first_char){
				conserved = false;
				break;
			}
		}
		if(conserved){
			out << "*";
		}
		else{
			out << " ";
		}
	}
	out << endl;

	return out.str();
}

string Alignment::ClustalWFormat(const vector<int> & columns){
	stringstream out;
	out << "CLUSTAL W\n\n";

	for(unsigned int i = 0; i < alignedSeqs.size();i++){
		stringstream temp;
		temp << i;
		assert(temp.str().length() <= 10);
		//names are up to 10 chars long, followed by a space
		for(unsigned int j = temp.str().length(); j < 11; j++){
			temp << " ";
		}

		//print only columns that we are interested in
		for(unsigned int j = 0; j < columns.size(); j++){
			temp << alignedSeqs[i]->at(columns[j]);
		}
		//temp << *alignedSeqs[i];
		out << temp.str() << endl;
	}
	//print line with *s
	for(unsigned int j = 0; j < 11; j++){
		out << " ";
	}
	for(unsigned int j = 0; j < alignedStruct.size(); j++){
		bool conserved = true;
		char first_char = alignedSeqs[0]->at(j);
		for(unsigned int i = 1; i < alignedSeqs.size();i++){
			if(alignedSeqs[i]->at(j) != first_char){
				conserved = false;
				break;
			}
		}
		if(conserved){
			out << "*";
		}
		else{
			out << " ";
		}
	}
	out << endl;

	return out.str();
}

void Alignment::printHeatMap(Tree & tree){

	cout << "Log-likelihood ratio matrix\n";
	for(int i = 0; i < (int)alignedStruct.size(); i++){

		for (int j = 0; j < (int)alignedStruct.size(); j++){
			if(j != 0){
				cout << "\t";
			}
			if (j > i){
				cout << logLikelihood(i,j, tree);
			}
			else{
				cout << "0";
			}
		}
		cout << endl;

	}
}

void Alignment::sparseBpTable(Tree & tree){


	double treeLength = tree.totalLength();
	for(unsigned int i = 0; i < alignedStruct.size()-MIN_DIST; i++){
		for(int j = i + MIN_DIST; j < (int)alignedStruct.size(); j++){
			int featureCount = 1;
			if(alignedStruct[i] == j){
				cout << "+1 ";
			}
			else{
				cout << "-1 ";
			}
			cout << featureCount++ << ":" << logPairedLikelihood(i,j, tree) / treeLength << " ";
			cout << featureCount++ << ":" << (logUnpairedLikelihood(i, tree) + logUnpairedLikelihood(j, tree)) / treeLength << " ";

			cout << featureCount++ << ":" << logPairedLikelihood(i+1,j-1, tree) / treeLength << " ";
			cout << featureCount++ << ":" << (logUnpairedLikelihood(i+1, tree) + logUnpairedLikelihood(j-1, tree)) / treeLength << " ";

			if(i > 0 && j+1 < (int)alignedStruct.size()){
				cout << featureCount++ << ":" << logPairedLikelihood(i-1,j+1, tree) / treeLength << " ";
				cout << featureCount++ << ":" << (logUnpairedLikelihood(i-1, tree) + logUnpairedLikelihood(j+1, tree)) / treeLength << " ";
			}
			else{
				featureCount += 2;
			}
			cout << featureCount++ << ":" << j-i << " ";
			cout << featureCount++ << ":" << consensusBPPercent(i, j) << " ";
			cout << featureCount++ << ":" << getCovariance(i,j);

			cout << endl;
		}
	}
}

string Alignment::FastaFormat(bool aligned) const{
	string out = "";
	for(unsigned int i = 0; i < alignedSeqs.size(); i++){
		out += ">" + *seqNames[i] + "\n";
		if(aligned){
			out += *alignedSeqs[i] + "\n";
		}
		else{
			out += *seqs[i] + "\n";
		}
	}
	return out;
}

void Alignment::columnTable(Tree & tree){
	cerr << "not yet implemented...\n";
	exit(-1);
}

void Alignment::emptyStruct(){
	alignedStruct = vector<int>(alignedSeqs[0]->length(), -1);
	fillSeqStructs();
}

