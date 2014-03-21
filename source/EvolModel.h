/*
 * Copyright: Nicholas P Wiebe (2008-2009) and Irmtraud M Meyer (2008-2009)
 * License: licensed under the GNU General Public License version 3 (GPLv3)
 * EvolModel.h
 * Singleton object that stores the models sequence evolution in paired and unpaired RNA
 * (Same model as in SimulFold, PFold)
 */

#ifndef EVOLMODEL_H_
#define EVOLMODEL_H_


class EvolModel {
public:
	virtual ~EvolModel();


	static double ePiSingle[4];
	static double eVSingle[4][4];
	static double eDSingle[4];
	static double eWSingle[4][4];

	static double ePiDouble[16];
	static double eVDouble[16][16];
	static double eDDouble[16];
	static double eWDouble[16][16];



private:
	EvolModel(); //singleton
};

#endif /* EVOLMODEL_H_ */
