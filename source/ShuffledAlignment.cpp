/*
 * Copyright: Nicholas P Wiebe (2008-2009) and Irmtraud M Meyer (2008-2009)
 * License: licensed under the GNU General Public License version 3 (GPLv3)
 */

#include "ShuffledAlignment.h"
#include <cassert>
#include <list>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include "Utilities.h"

//const string ShuffledAlignment::TEMP_FASTA_FILENAME = "temp.fasta";
//const string ShuffledAlignment::TEMP_OUTPUT_FASTA_FILENAME = "temp.out.fasta";
//const string ShuffledAlignment::REALIGNER_LOCATION = "Realigner/bin";

//const string ShuffledAlignment::RNAZ_SHUFFLER_LOC = "RNAz_perl/rnazRandomizeAln.pl";

ShuffledAlignment::ShuffledAlignment(const Alignment & a){

	alignmentName = a.alignmentName;//TODO: get rid of alignment name concept

	for(unsigned int i = 0; i < a.seqNames.size(); i++){
		seqNames.push_back(new string(*a.seqNames[i]));
	}



	//copy alignment
	for(unsigned int i = 0; i < a.alignedSeqs.size(); i++){
		alignedSeqs.push_back(new string(*(a.alignedSeqs[i])));
	}

	//clear structure
	for(unsigned int i = 0; i < a.alignedStruct.size();i++){
		alignedStruct.push_back(-1);
	}

	//use RNAz to shuffle columns
	shuffleRNAz();


//
//	vector<unsigned int> shuffleableColumns;
//
//	//copy struct, preserving only bps in 'th'
//	//find columns that are not part of helix 'th':
//	assert(a.alignmentLabels.size() == a.alignedStruct.size());
//	for(unsigned int i = 0; i < a.alignedStruct.size(); i++){
//		if(a.alignedStruct[i] != -1){
//			alignedStruct.push_back(a.alignedStruct[i]);
//		}
//		else{
//			alignedStruct.push_back(-1);
//			shuffleableColumns.push_back(i);
//		}
//	}
//
//
//	//cout << "number of unpaired columns\n";
//	//cout << "The value of RAND_MAX is " <<  RAND_MAX << endl;
//
//	//shuffle list unpaired columns:
//	for(unsigned int i = 0; i < shuffleableColumns.size(); i++){
//		unsigned int random_column = rand();
//
//		random_column = i + random_column % (shuffleableColumns.size()-i);
//
//		//swap i with random_column
//		int temp = shuffleableColumns[i];
//		shuffleableColumns[i] = shuffleableColumns[random_column];
//		shuffleableColumns[random_column] = temp;
//	}
//
//	//reorder columns:
//
//
//	int shuffleableColsIndex = 0;
//	for(unsigned int i = 0; i < alignedStruct.size(); i++){
//		if(alignedStruct[i] == -1){
//			for(unsigned int j = 0; j < a.alignedSeqs.size(); j++){
//				alignedSeqs[j]->at(i) = a.alignedSeqs[j]->at(shuffleableColumns[shuffleableColsIndex]);
//			}
//			shuffleableColsIndex++;
//		}
//	}
//

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
	labelHelices2();
	initializeStatsMatrix();

}

ShuffledAlignment::ShuffledAlignment(const Alignment & a, int th, string treeFile)
{
	alignmentName = a.alignmentName;//TODO: get rid of alignment name concept

	for(unsigned int i = 0; i < a.seqNames.size(); i++){
		seqNames.push_back(new string(*a.seqNames[i]));
	}


	vector<string> tempAlignment;

	//sanity check:
	assert(a.alignmentLabels.size() == a.alignedStruct.size());

	//if realign is true, produce new alignment blocks and
	//update struct vector
	if(treeFile.compare("") != 0){

		//initialize new alignment strings
		unsigned int i;
		for(i = 0; i < a.alignedSeqs.size(); i++){
			tempAlignment.push_back(string(""));
		}

		vector<int> offset_vec;
		int offset = 0;
		int realign_start = 0;

		for(i = 0; i < a.alignedStruct.size(); i++){
			if(a.alignmentLabels[i] == th){
				if((int)i > realign_start){

					vector<string> interval = a.realignInterval(realign_start, i, treeFile);

					//add new alignment block
					for(unsigned int j = 0; j < interval.size(); j++){
						tempAlignment[j] += interval[j];
					}
					for(unsigned int j = 0; j < interval[0].length(); j++){
						alignedStruct.push_back(-1);
					}
					offset += interval[0].length() - (i - realign_start);

				}
				//add current column to alignment
				for(unsigned int j = 0; j < tempAlignment.size(); j++){
					tempAlignment[j] += a.alignedSeqs[j]->substr(i, 1);
				}

				offset_vec.push_back(offset);
				if(a.alignedStruct[i] < (int)i){
					assert(alignedStruct[a.alignedStruct[i] + offset_vec[a.alignedStruct[i]]] == -2);
					alignedStruct.push_back(a.alignedStruct[i] + offset_vec[a.alignedStruct[i]]);
					alignedStruct[a.alignedStruct[i] + offset_vec[a.alignedStruct[i]]] = alignedStruct.size() - 1;
				}
				else{
					alignedStruct.push_back(-2); //will be set when 3' pairing partner found
				}
				realign_start = i+1;
			}
			else{
				offset_vec.push_back(offset);
			}
		}
		//get last alignment block (if needed)
		if((int)i > realign_start){
			vector<string> interval = a.realignInterval(realign_start, i, treeFile);

			//add new alignment block
			for(unsigned int j = 0; j < interval.size(); j++){
				tempAlignment[j] += interval[j];
			}
			for(unsigned int j = 0; j < interval[0].length(); j++){
				alignedStruct.push_back(-1);
			}
			offset += interval[0].length() - (i - realign_start);

		}

		//cout << tempAlignment[0] << endl << *a.alignedSeqs[0] << endl;

		//sanity checks:
		for(i = 0; i < tempAlignment.size(); i++){
			assert(tempAlignment[i].length() == tempAlignment[0].length());
		}

		for(i = 0; i < alignedStruct.size(); i++){
			assert(alignedStruct[i] >= -1);
			if(alignedStruct[i] >=0){
				assert((int)i == alignedStruct[alignedStruct[i]]);
			}
		}

	}
	else{
		//copy alignment
		for(unsigned int i = 0; i < a.alignedSeqs.size(); i++){
			tempAlignment.push_back(*(a.alignedSeqs[i]));
		}

		//copy struct, preserving only bps in 'th'
		for(unsigned int i = 0; i < a.alignedStruct.size(); i++){
			if(a.alignmentLabels[i] == th){
					alignedStruct.push_back(a.alignedStruct[i]);
			}
			else{
				alignedStruct.push_back(-1);
			}
		}

	}

	//find shuffleable columns
	vector<int> shuffleableColumns;
	for(unsigned int i = 0; i < alignedStruct.size(); i++){

		if(alignedStruct[i] == -1){
			shuffleableColumns.push_back(i);
		}
	}

	//initialize alignedSeqs
	for(unsigned int i = 0; i < a.alignedSeqs.size(); i++){
		alignedSeqs.push_back(new string(tempAlignment[i]));
	}

	shuffleRNAz(shuffleableColumns);
//
//	//cout << "number of unpaired columns\n";
//
//	srand((unsigned)time(0));
//	//cout << "The value of RAND_MAX is " <<  RAND_MAX << endl;
//
//	//shuffle list unpaired columns:
//	for(unsigned int i = 0; i < shuffleableColumns.size(); i++){
//		unsigned int random_column = rand();
//
//		random_column = i + random_column % (shuffleableColumns.size()-i);
//
//		//swap i with random_column
//		int temp = shuffleableColumns[i];
//		shuffleableColumns[i] = shuffleableColumns[random_column];
//		shuffleableColumns[random_column] = temp;
//	}
//
//	//reorder columns:
//
//	for(unsigned int i = 0; i < a.alignedSeqs.size(); i++){
//		alignedSeqs.push_back(new string(""));
//	}
//
//	int shuffleableColsIndex = 0;
//	for(unsigned int i = 0; i < alignedStruct.size(); i++){
//		if(alignedStruct[i] == -1){
//			for(unsigned int j = 0; j < tempAlignment.size(); j++){
//				*alignedSeqs[j] += tempAlignment[j].at(shuffleableColumns[shuffleableColsIndex]);
//			}
//			shuffleableColsIndex++;
//		}
//		else{
//			for(unsigned int j = 0; j < tempAlignment.size(); j++){
//				*alignedSeqs[j] += tempAlignment[j].at(i);
//			}
//		}
//	}


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
	labelHelices2();
	initializeStatsMatrix();

	//sanity check:
	if(treeFile.compare("") != 0){
		for(unsigned int i = 0; i < seqs.size(); i++){
			//cout << i << endl;
			//cout << seqs[i]->length() << "\t" << a.seqs[i]->length() << endl;
			//cout << *seqs[i] << endl << *a.seqs[i] << endl;

			assert(seqs[i]->length() == a.seqs[i]->length());
			assert(helixLabels[i]->size() == a.helixLabels[i]->size());


			unsigned int j = 0;
			unsigned int k = 0;
			while(true){
				while(j < helixLabels[i]->size() && helixLabels[i]->at(j) != 0){
					j++;
				}
				while(k < a.helixLabels[i]->size() && a.helixLabels[i]->at(k) != th){
					k++;
				}
				if(k < a.helixLabels[i]->size() && j < helixLabels[i]->size()){
					assert(a.seqs[i]->at(k) == seqs[i]->at(j));
				}
				else{
					assert(k >= a.helixLabels[i]->size() && j >= helixLabels[i]->size());
					//cout << "getting out!\n";
					break;
				}
				j++;
				k++;
			}
		}
	}

}

ShuffledAlignment::~ShuffledAlignment() {
	clearAll();
}

void ShuffledAlignment::shuffleRNAz(){

	unsigned int original_length = alignedSeqs[0]->length();

	//write clustal file:
	string clustalFilenameTemplate = "/tmp/clustal.XXXXXX";
	char clustalFilename[clustalFilenameTemplate.size()+1];

	strcpy(clustalFilename, clustalFilenameTemplate.c_str());

	int unique_fd;
	FILE * unique_file;

	if ((unique_fd = mkstemp(clustalFilename)) == -1 ){
		cerr << "Could not create temporary clustal file for shuffling\n";
		exit(-1);
	}

	if((unique_file = fdopen(unique_fd, "w")) == NULL){
		cerr << "Could not open temporary clustal file for shuffling\n";
		exit(-1);
	}

	fputs(ClustalWFormat().c_str(),unique_file);
	fclose(unique_file);

	//make a temp file for output:
	string outTemplate = "/tmp/shuffle_out.XXXXXX";
	char outFilename[outTemplate.size()+1];
	strcpy(outFilename, outTemplate.c_str());

	if ((unique_fd = mkstemp(outFilename)) == -1 ){
		cerr << "Could not create temporary out file for shuffling\n";
		exit(-1);
	}

	if((unique_file = fdopen(unique_fd, "w")) == NULL){
		cerr << "Could not open temporary out file for shuffling\n";
		exit(-1);
	}
	string temp_contents = "TEMP!\n"; // this will be replaced by rnazShuffler output
	fputs(temp_contents.c_str(), unique_file);
	fclose(unique_file);

	//run shuffler:
	string command = string(RNAZ_SHUFFLER_LOC) + " -l 1 " + string(clustalFilename) + " > " + string(outFilename);
	if(system(command.c_str()) != 0){
		cerr << "Error: RNAz Shuffler encountered an error\n";
		exit(-1);
	}

	//read output:
	parseClustalWFile(string(outFilename));

	assert(original_length == alignedSeqs[0]->length());

	//clean up files
	command = "rm -f " + string(outFilename);
	if(system(command.c_str()) != 0){
		cerr << "Error cleaning up temporary shuffle output file\n";
		exit(-1);
	}

	command = "rm -f " + string(clustalFilename);
	if(system(command.c_str()) != 0){
		cerr << "Error cleaning up temporary shuffle input file\n";
		exit(-1);
	}

}

void ShuffledAlignment::parseClustalWFile(const string & filename){

	//clear current alignment:
	for(unsigned int i = 0; i < alignedSeqs.size(); i++){
		*alignedSeqs[i] = "";
	}

	ifstream shuffleFile(filename.c_str());
	if(!shuffleFile.is_open()){
		cerr << "Error: cannot open clustal file " << filename << endl;
		exit(-1);
	}
	string line;
	getline(shuffleFile, line);

	if(line.substr(0, 9).compare("CLUSTAL W")){
		cerr << "Error: file missing CLUSTAL W header\n";
	}

	while(shuffleFile.good()){
		getline(shuffleFile, line);
		Utilities::TrimSpaces(line);
		if(line.compare("")){
			//start of alignment block:
			int start;
			for(unsigned int i = 0; i < alignedSeqs.size() - 1; i++){
				start = line.find_last_of(" \t");
				*alignedSeqs[i] = *alignedSeqs[i] + line.substr(start+1);

				assert(shuffleFile.good());
				getline(shuffleFile,line);
				Utilities::TrimSpaces(line);
				assert(line.compare("")); //not blank
			}
			//last sequence in a block
			start = line.find_last_of(" \t");
			*alignedSeqs[alignedSeqs.size() - 1] = *alignedSeqs[alignedSeqs.size() - 1] + line.substr(start+1);
			//make sure we're at the end of a block
			if(shuffleFile.good()){
				getline(shuffleFile,line);
				Utilities::TrimSpaces(line);
				assert(!line.compare("")); //blank
			}
		}
	}

}

void ShuffledAlignment::shuffleRNAz(const vector<int> & columns){

	unsigned int original_length = alignedSeqs[0]->length();

	//write clustal file:
	string clustalFilenameTemplate = "/tmp/clustal.XXXXXX";
	char clustalFilename[clustalFilenameTemplate.size()+1];

	strcpy(clustalFilename, clustalFilenameTemplate.c_str());

	int unique_fd;
	FILE * unique_file;

	if ((unique_fd = mkstemp(clustalFilename)) == -1 ){
		cerr << "Could not create temporary clustal file for shuffling\n";
		exit(-1);
	}

	if((unique_file = fdopen(unique_fd, "w")) == NULL){
		cerr << "Could not open temporary clustal file for shuffling\n";
		exit(-1);
	}

	fputs(ClustalWFormat(columns).c_str(),unique_file);
	fclose(unique_file);

	//make a temp file for output:
	string outTemplate = "/tmp/shuffle_out.XXXXXX";
	char outFilename[outTemplate.size()+1];
	strcpy(outFilename, outTemplate.c_str());

	if ((unique_fd = mkstemp(outFilename)) == -1 ){
		cerr << "Could not create temporary out file for shuffling\n";
		exit(-1);
	}

	if((unique_file = fdopen(unique_fd, "w")) == NULL){
		cerr << "Could not open temporary out file for shuffling\n";
		exit(-1);
	}
	string temp_contents = "TEMP!\n"; // this will be replaced by rnazShuffler output
	fputs(temp_contents.c_str(), unique_file);
	fclose(unique_file);

	//run shuffler:
	string command = string(RNAZ_SHUFFLER_LOC) + " -l 1 " + string(clustalFilename) + " > " + string(outFilename);
	if(system(command.c_str()) != 0){
		cerr << "Error: RNAz Shuffler encountered an error\n";
		exit(-1);
	}

	//read output:
	parseClustalWFile(string(outFilename), columns);

	assert(original_length == alignedSeqs[0]->length());

	//clean up files
	command = "rm -f " + string(outFilename);
	if(system(command.c_str()) != 0){
		cerr << "Error cleaning up temporary shuffle output file\n";
		exit(-1);
	}

	command = "rm -f " + string(clustalFilename);
	if(system(command.c_str()) != 0){
		cerr << "Error cleaning up temporary shuffle input file\n";
		exit(-1);
	}

}

void ShuffledAlignment::parseClustalWFile(const string & filename, const vector<int> & columns){

	ifstream shuffleFile(filename.c_str());
	if(!shuffleFile.is_open()){
		cerr << "Error: cannot open clustal file " << filename << endl;
		exit(-1);
	}
	string line;
	getline(shuffleFile, line);

	if(line.substr(0, 9).compare("CLUSTAL W")){
		cerr << "Error: file missing CLUSTAL W header\n";
	}

	int blockOffset = 0;
	string seq;
	while(shuffleFile.good()){
		getline(shuffleFile, line);
		Utilities::TrimSpaces(line);
		if(line.compare("")){
			//start of alignment block:
			int start;
			for(unsigned int i = 0; i < alignedSeqs.size() - 1; i++){
				start = line.find_last_of(" \t");
				seq = line.substr(start+1);

				//copy columns into correct locations
				assert(blockOffset + seq.length() - 1 < columns.size());
				for(unsigned int j = 0; j < seq.length(); j++){
					alignedSeqs[i]->at(columns[blockOffset + j]) = seq[j];
				}

				assert(shuffleFile.good());
				getline(shuffleFile,line);
				Utilities::TrimSpaces(line);
				assert(line.compare("")); //not blank
			}
			//last sequence in a block
			start = line.find_last_of(" \t");

			seq = line.substr(start+1);
			assert(blockOffset + seq.length() - 1 < columns.size());
			for(unsigned int j = 0; j < seq.length(); j++){
				alignedSeqs[alignedSeqs.size()-1]->at(columns[blockOffset + j]) = seq[j];
			}

			blockOffset += seq.length();
			//make sure we're at the end of a block
			if(shuffleFile.good()){
				getline(shuffleFile,line);
				Utilities::TrimSpaces(line);
				assert(!line.compare("")); //blank
			}
		}
	}
}

