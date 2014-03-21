/*
 * Copyright: Nicholas P Wiebe (2008-2009) and Irmtraud M Meyer (2008-2009)
 * License: licensed under the GNU General Public License version 3 (GPLv3)
 * AlignedHelix.cpp
 */

#include "AlignedHelix.h"
#include <math.h>
#include "Utilities.h"
#include "Alignment.h"
#include <cassert>
#include <sstream>

AlignedHelix::AlignedHelix(const Helix & helix, int seqIndex, vector<int> & seq2AlignmentMap) {

	int alignedPos5, alignedPos3;

	for(int i = 0; i < helix.length; i++){
		alignedPos5 = seq2AlignmentMap[helix.pos5+i];
		alignedPos3 = seq2AlignmentMap[helix.pos3-i];

		bps.push_back(pair<int,int>(alignedPos5, alignedPos3));
	}

	appearsIn.insert(seqIndex);
}

AlignedHelix::AlignedHelix(int outerBpPos5, int outerBpPos3, int length){
	assert(outerBpPos5 >= 0);

	int pos5, pos3;
	for(int i = 0; i < length; i++){
		pos5 = outerBpPos5 + i;
		pos3 = outerBpPos3 - i;
		assert(pos3 > pos5);
		bps.push_back(pair<int,int>(pos5, pos3));
	}
}


AlignedHelix::~AlignedHelix() {

}

string AlignedHelix::dotBracket(unsigned int alignmentLength){
	string out = "";
	for(unsigned int i = 0; i < alignmentLength; i++){
		out += ".";
	}

	for(unsigned int i = 0; i < bps.size(); i++){
		out[bps[i].first] = '(';
		out[bps[i].second] = ')';
	}

	return out;
}

bool AlignedHelix::insert(const Helix & helix, int seqIndex, vector<int> & seq2AlignmentMap){
	int alignedPos5, alignedPos3;

	//has to match exactly, including size!
	if(helix.length != (int)bps.size()){
		return false;
	}

	for(int i = 0; i < helix.length; i++){
		alignedPos5 = seq2AlignmentMap[helix.pos5+i];
		alignedPos3 = seq2AlignmentMap[helix.pos3-i];

		if(alignedPos5 != bps[i].first || alignedPos3 != bps[i].second){
			return false;
		}
	}
	appearsIn.insert(seqIndex);
	return true;
}

double AlignedHelix::logLikeRatio(Alignment & a, Tree & tree){
	int pos3, pos5;
	double loglike = 0;
	for(vector<pair<int, int> >::iterator it = bps.begin(); it != bps.end(); it++){
		pos5 = it->first;
		pos3 = it->second;
		loglike += a.logLikelihood(pos5,pos3, tree);
	}
	return loglike / bps.size();
}

double AlignedHelix::pairedLikelihood(Alignment & a, Tree & tree){
	int pos3, pos5;
	double paired = 0;
	for(vector<pair<int, int> >::iterator it = bps.begin(); it != bps.end(); it++){
		pos5 = it->first;
		pos3 = it->second;
		paired += log(tree.calcFelsDouble(a, pos5, pos3))/log(2);
	}
	return paired / bps.size();
}
double AlignedHelix::unpairdLikelihood(Alignment & a, Tree & tree){
	int pos3, pos5;
	double unpaired = 0;
	for(vector<pair<int, int> >::iterator it = bps.begin(); it != bps.end(); it++){
		pos5 = it->first;
		pos3 = it->second;
		unpaired += log(tree.calcFelsSingle(a, pos5))/log(2);
		unpaired += log(tree.calcFelsSingle(a, pos3))/log(2);
	}
	return unpaired / bps.size();
}

double AlignedHelix::canonicalBP(Alignment & a){
	int pos3, pos5;
	double canonicalBP = 0;

	for(vector<pair<int, int> >::iterator it = bps.begin(); it != bps.end(); it++){
		pos5 = it->first;
		pos3 = it->second;

		canonicalBP += a.consensusBPPercent(pos5, pos3);
	}
	return canonicalBP / bps.size();
}

double AlignedHelix::covariance(Alignment & a){
	int pos3, pos5;
	double cov = 0;

	for(vector<pair<int, int> >::iterator it = bps.begin(); it != bps.end(); it++){
		pos5 = it->first;
		pos3 = it->second;

		cov += a.getCovariance(pos5, pos3);
	}
	return cov / bps.size();
}

double AlignedHelix::conservation(Alignment & a){
	int pos3, pos5;
	double conserve = 0;

	for(vector<pair<int, int> >::iterator it = bps.begin(); it != bps.end(); it++){
		pos5 = it->first;
		pos3 = it->second;

		conserve += a.getSeqCons(pos5);
		conserve += a.getSeqCons(pos3);
	}
	return conserve / (2*bps.size());
}

bool AlignedHelix::isConsensusHelix(Alignment & a){
	int pos3,pos5;

	for(unsigned int i = 0; i < bps.size(); i++){

		pos5 = bps[i].first;
		pos3 = bps[i].second;

		if(a.alignedStruct[pos5] != pos3){
			return false;
		}

		if(a.alignedStruct[pos5 + 1] == pos3-1){
			//bp inner to current bp is in consensus structure
			if(!(i + 1 < bps.size()) || !(bps[i+1] == pair<int, int>(pos5+1, pos3-1))){
				//bp inner to current bp not in this helix
				return false;
			}
		}

		if((pos5 > 0 && pos3 + 1 < (int)a.alignedStruct.size()) && (a.alignedStruct[pos5 - 1] == pos3+1)){
			//bp outer to current bp is in consensus structure
			if(!(i > 0) || !(bps[i-1] == pair<int, int>(pos5-1, pos3+1))){
				return false;
			}
		}
	}

	return true;
}

bool AlignedHelix::isPartialConsensusHelix(Alignment & a){
	int pos3, pos5;

	for(unsigned int i = 0; i < bps.size(); i++){

		pos5 = bps[i].first;
		pos3 = bps[i].second;

		if(a.alignedStruct[pos5] == pos3){
			return true;
		}
	}

	return false;
}

bool AlignedHelix::isCompetingHelix(Alignment & a){
	int pos3, pos5;

	for(unsigned int i = 0; i < bps.size(); i++){

		for(set<int>::iterator seqIt = appearsIn.begin(); seqIt != appearsIn.end(); seqIt++){
			pos5 = a.alignment2SeqMap[*seqIt]->at(bps[i].first);
			pos3 = a.alignment2SeqMap[*seqIt]->at(bps[i].second);

			//sanity check: make sure that positions are not gaps in the sequence
			assert(pos5>=0);
			assert(pos3>=0);

			//note: maybe I should check if helix is partial true helix?

			if(a.seqStructs[*seqIt]->at(pos5) > 0 && a.seqStructs[*seqIt]->at(pos5) != pos3){
				//5' position competing
				if(Utilities::iCheckLongHelix(*a.seqStructs[*seqIt], pos5)){
					return true;
				}
			}
			if(a.seqStructs[*seqIt]->at(pos3) > 0 && a.seqStructs[*seqIt]->at(pos3) != pos5){
				//3' position competing
				if(Utilities::iCheckLongHelix(*a.seqStructs[*seqIt], pos3)){
					return true;
				}
			}

		}

	}

	return false;
}

int AlignedHelix::consensusBps(Alignment & a){
	int pos3, pos5;
	int consensusBpsCount = 0;
	for(unsigned int i = 0; i < bps.size(); i++){

		pos5 = bps[i].first;
		pos3 = bps[i].second;

		if(a.alignedStruct[pos5] == pos3){
			consensusBpsCount++;
		}
	}

	return consensusBpsCount;
}

double AlignedHelix::midpoint(){
	int positionTotal = 0;
	for(unsigned int i = 0; i < bps.size(); i++){
		positionTotal += bps[i].first;
		positionTotal += bps[i].second;
	}

	return positionTotal / ((double) bps.size() * 2);
}

double AlignedHelix::midpoint(Alignment & a){
	return midpoint() / a.alignedStruct.size();
}

int AlignedHelix::length(){
	return bps.size();
}


int AlignedHelix::appearances(){
	return appearsIn.size();
}

string AlignedHelix::bpsString(){
	stringstream out;
	unsigned int i;
	for(i = 0; i < bps.size() - 1; i++){
		out << bps[i].first << ":" << bps[i].second << ",";
	}
	out << bps[i].first << ":" << bps[i].second;
	return out.str();
}


pair<double,double> AlignedHelix::cisTransScore(Alignment & a){
	pair<double, double> score(0,0);

	//Cis/Trans Definitions
	//c pairs with i in the competing helix
	//i pairs with i^ in the consensus structure
	//c:i:i^ = 5'-cis
	//c:i^:i = 5'-trans
	//i^:i:c = 3'-cis
	//i:i^:c = 3'-trans

	int i_bar, i, c;
	double temp;
	int L = a.alignedStruct.size();

	for (unsigned int j = 0; j < bps.size(); j++){

		if(a.alignedStruct[bps[j].first] >= 0){
			i_bar = a.alignedStruct[bps[j].first];
			c = bps[j].second;
			i = bps[j].first;
			if(c != i_bar){
				if(i_bar < i){
					//3'-cis
					temp = 1/((c-i) * log(L-i));
					assert(temp > 0);
					score.first -= temp;
					//score.first -= 1/((c-i) * log(L-i));
				}
				else if(i_bar < c){
					//i_bar > i
					//3'-trans
					temp = 1/((c-i_bar) * log(L-i_bar));
					assert(temp > 0);
					score.second += temp;
					//score.second += 1/((c-i_bar) * log(L-i_bar));
				}
				//else 3'-mid : ignore
			}
		}

		if(a.alignedStruct[bps[j].second] >= 0){
			i_bar = a.alignedStruct[bps[j].second];
			c = bps[j].first;
			i = bps[j].second;
			if(c != i_bar){
				if(i_bar > i){
					//5'-cis
					temp = 1/((i-c) * log(i+1));
					assert(temp > 0);
					score.first += temp;
					//score.first += 1/((i-c) * log(i));
				}
				else if(c < i_bar){
					//i_bar < i
					//5'-trans
					temp =1/((i_bar - c) * log(i_bar+1));
					assert(temp > 0);
					score.second -= temp;
					//score.second -= 1/((i_bar - c) * log(i_bar));
				}
				//need the +1 in log term to make it symmetrical
				//else 5'mid : ignore
			}

		}

	}

	return score;
}


void AlignedHelix::competeScore(Alignment & a, double & cis5, double & cis3, double & trans5, double & trans3, double & mid5, double & mid3){

	cis5 = 0;
	cis3 = 0;
	trans5 = 0;
	trans3 = 0;
	mid5 = 0;
	mid3 = 0;

	//Cis/Trans Definitions
	//c pairs with i in the competing helix
	//i pairs with i^ in the consensus structure
	//c:i:i^ = 5'-cis
	//c:i^:i = 5'-trans
	//i^:i:c = 3'-cis
	//i:i^:c = 3'-trans
	//i^:c:i = 5' mid --> competing base is between i and i^, and 5' of the base it pairs with in the competing helix
	//i:c:i^ = 3' mid --> competing base is between i and i^, and 3' of the base it pairs with in the competing helix


	int i_bar, i, c;
	double temp;
	int L = a.alignedStruct.size();

	for (unsigned int j = 0; j < bps.size(); j++){

		if(a.alignedStruct[bps[j].first] >= 0){
			i_bar = a.alignedStruct[bps[j].first];
			c = bps[j].second;
			i = bps[j].first;
			if(c != i_bar){
				if(i_bar < i){
					//3'-cis
					temp = 1/((c-i) * log(L-i));
					assert(temp > 0);
					cis3 += temp;
					//score.first -= 1/((c-i) * log(L-i));
				}
				else if(i_bar < c){
					//i_bar > i
					//3'-trans
					temp = 1/((c-i_bar) * log(L-i_bar));
					assert(temp > 0);
					trans3 += temp;
					//score.second += 1/((c-i_bar) * log(L-i_bar));
				}
				else{
					//else 3'-mid
					temp = 1 / ((c-i) * log(i_bar - i));
					assert(temp > 0);
					mid3 += temp;
				}

			}
		}

		if(a.alignedStruct[bps[j].second] >= 0){
			i_bar = a.alignedStruct[bps[j].second];
			c = bps[j].first;
			i = bps[j].second;
			if(c != i_bar){
				if(i_bar > i){
					//5'-cis
					temp = 1/((i-c) * log(i+1));
					assert(temp > 0);
					cis5 += temp;
					//score.first += 1/((i-c) * log(i));
				}
				else if(c < i_bar){
					//i_bar < i
					//5'-trans
					temp =1/((i_bar - c) * log(i_bar+1));
					assert(temp > 0);
					trans5 += temp;
					//score.second -= 1/((i_bar - c) * log(i_bar));
				}
				//need the +1 in log term to make it symmetrical
				else{
					//else 5'mid
					temp = 1 / ((i-c) * log(i - i_bar));
					assert(temp > 0);
					mid5 += temp;
				}


			}

		}

	}
}
