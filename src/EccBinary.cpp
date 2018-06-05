/*
 * EccBinary.cpp
 *
 *  Created on: Jun 5, 2018
 *      Author: Travis Robson
 */

#include "EccBinary.h"

//namespace eb {

EccBinary::EccBinary(double m1, double m2, double e0, double F0) {
	// TODO Auto-generated constructor stub
	this->m1 = m1;
	this->m2 = m2;
	this->e0 = e0;
	this->F0 = F0;

	//this->mu = calc_mu();
}

EccBinary::~EccBinary() {
	// TODO Auto-generated destructor stub
}

//void EccBinary::calc_mu()
//{
//
//	return;
//}

//} /* namespace eb */
