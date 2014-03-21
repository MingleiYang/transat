/*
 * Copyright: Nicholas P Wiebe (2008-2009) and Irmtraud M Meyer (2008-2009)
 * License: licensed under the GNU General Public License version 3 (GPLv3)
 */

#include "StatsWrapper.h"
#include <iostream>

StatsWrapper::StatsWrapper() {

	iStats = new int[3];
	revIStats = new int[3];

	for(int i = 0; i < 3; i++){
		iStats[i] = 0;
		revIStats[i] = 0;
	}

	eStats = new double[3];
	revEStats = new double[3];

	for(int i = 0; i < 3; i++){
		eStats[i] = 0.0;
		revEStats[i] = 0.0;
	}

	pStat = 0.0;
	gStat = 0.0;
	revgStat = 0.0;
	type = NONE;
	revtype = NONE;
}

StatsWrapper::StatsWrapper(const StatsWrapper & other)
{
	iStats = new int[3];
	revIStats = new int[3];

	for(int i = 0; i < 3; i++){
		iStats[i] = other.iStats[i];
		revIStats[i] = other.revIStats[i];
	}

	eStats = new double[3];
	revEStats = new double[3];

	for(int i = 0; i < 3; i++){
		eStats[i] = other.eStats[i];
		revEStats[i] = other.revEStats[i];

	}

	pStat = other.pStat;
	gStat = other.gStat;
	revgStat = other.revgStat;
	type = other.type;
	revtype = other.revtype;
}

StatsWrapper::~StatsWrapper() {

	//std::cout << "Destructor!\n";
	delete[] eStats;
	delete[] revEStats;
	delete[] iStats;
	delete[] revIStats;
}

StatsWrapper & StatsWrapper::operator+=(const StatsWrapper &rhs){
	//sum up eStats and iStats
	//don't bother with pStat, gStat, etc.

	for(int i = 0; i < 3; i++){
		iStats[i] += rhs.iStats[i];
		revIStats[i] += rhs.revIStats[i];
		eStats[i] += rhs.eStats[i];
		revEStats[i] += rhs.revEStats[i];
	}

	return *this;
}

double StatsWrapper::mid3(){
	return eStats[2];
}

double StatsWrapper::mid5(){
	return revEStats[2];
}
double StatsWrapper::cis5()
{
	return revEStats[0];
}
double StatsWrapper::trans5()
{
	return revEStats[1];
}
double StatsWrapper::cis3()
{
	return eStats[0];
}
double StatsWrapper::trans3()
{
	return eStats[1];
}
double StatsWrapper::cis()
{
	return cis5() - cis3();
}
double StatsWrapper::trans()
{
	return trans3() - trans5();
}
double StatsWrapper::mid(){
	return mid5() - mid3();
}

bool StatsWrapper::isZero() {
	if(
			( cis5()   == 0) && // cis5
			( cis3()   == 0) && // cis3
			( trans5() == 0) && // trans5
			( trans3() == 0) && // trans3
			( mid3() == 0) &&
			(mid5() == 0)
		) {
		return 1;
	} else {
		return 0;
	}
}

const StatsWrapper StatsWrapper::operator+(const StatsWrapper & other) const
{
	return StatsWrapper(*this) += other;
}


