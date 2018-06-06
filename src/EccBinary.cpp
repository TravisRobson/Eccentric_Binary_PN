/*
 * EccBinary.cpp
 *
 *  Created on: Jun 5, 2018
 *      Author: Travis Robson
 */

#include <cmath>


#include "EccBinary.h"

//namespace eb {

EccBinary::EccBinary(double m1, double m2, double e0, double F0) {
	// TODO Auto-generated constructor stub
	this->m1 = m1;
	this->m2 = m2;
	this->e0 = e0;
	this->F0 = F0;

	calc_mass_parameters();
}

EccBinary::~EccBinary() {
	// TODO Auto-generated destructor stub
}

void EccBinary::calc_m()
{
	this->m = this->m1 + this->m2;
}

void EccBinary::calc_mu()
{
	this->mu = this->m1*this->m2/(this->m1+this->m2);
}

void EccBinary::calc_eta()
{
	this->eta = this->m1*this->m2/(this->m1+this->m2)/(this->m1+this->m2);
}

void EccBinary::calc_Mc()
{
	this->Mc = pow(this->m1*this->m2, 3./5.)/pow(this->m1 + this->m2, 1./5.);
}

void EccBinary::calc_mass_parameters()
{
	calc_m();
	calc_mu();
	calc_eta();
	calc_Mc();
}

//} /* namespace eb */
