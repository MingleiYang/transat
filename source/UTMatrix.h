/*
 * Copyright: Nicholas P Wiebe (2008-2009) and Irmtraud M Meyer (2008-2009)
 * License: licensed under the GNU General Public License version 3 (GPLv3)
 * UTMatrix.h
 *
 * Data structure for efficient storage of upper triangular matrices
 * (re-written to be a bit less weird... njwiebe, Sept 2009)
 */

#ifndef UTMATRIX_H_
#define UTMATRIX_H_

class UTMatrix {
public:
	UTMatrix(unsigned int size);
	virtual ~UTMatrix();
	double at(unsigned int i, unsigned int j);
	void set(unsigned int i, unsigned int j, double val);

protected:
//	double * data;
	double ** data;
	unsigned int size;
};

#endif /* UTMATRIX_H_ */
