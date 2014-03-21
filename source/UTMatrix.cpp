/*
 * Copyright: Nicholas P Wiebe (2008-2009) and Irmtraud M Meyer (2008-2009)
 * License: licensed under the GNU General Public License version 3 (GPLv3)
 */

#include "UTMatrix.h"
#include <cassert>

UTMatrix::UTMatrix(unsigned int size_) {
	assert(size_ > 0);

//	data = new double[((size_ + 1)* size_) / 2];
	size = size_;
	data = new double*[size_];
	for(unsigned int i = 0; i < size_; i++){
		data[i] = new double[size_ - i];
	}


}


UTMatrix::~UTMatrix() {

	for(unsigned int i = 0; i < size; i++){
		delete[] data[i];
	}
	delete[] data;
}

double UTMatrix::at(unsigned int i, unsigned int j){
	assert(j < size);
	assert(i<=j);
//	return data[i * size + j - ((i+1)*i)/2];

	return data[i][j-i];
}

void UTMatrix::set(unsigned int i, unsigned int j, double val){
	assert(j < size);
	assert(i <=j);

//	data[i * size + j - ((i+1)*i)/2] = val;
	data[i][j-i] = val;
}
