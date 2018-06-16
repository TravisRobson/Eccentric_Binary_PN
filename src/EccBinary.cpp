/*
 * EccBinary.cpp
 *
 *  Created on: Jun 5, 2018
 *      Author: Travis Robson
 */

#include <cmath>
#include <string>
#include <fstream>


#include "EccBinary.h"
#include "Constants.h"

using namespace std;


EccBinary::EccBinary(double m1, double m2, double et0, double F0, int PN) {
	char buf[64];
	char tag[256];

	// TODO Auto-generated constructor stub
	this->m1  = m1;
	this->m2  = m2;
	this->et0 = et0;
	this->F0  = F0;
	this->PN  = PN;

	calc_mass_parameters();

	this->orb = new Orbit();
	this->orb->set_N(10);    // dummy for now to make sure I coded things properly

	// construct a unique identifying tag
	sprintf(tag, "%.1f", this->m1/TSUN);
	strcat(tag, "_");

	sprintf(buf, "%.1f", this->m2/TSUN);
	strcat(tag, buf);
	strcat(tag, "_");

	sprintf(buf, "%.2f", this->et0);
	strcat(tag, buf);
	strcat(tag, "_");

	sprintf(buf, "%.1f", this->F0);
	strcat(tag, buf);
	strcat(tag, "_");

	sprintf(buf, "%d", this->PN);
	strcat(tag, buf);
	strcat(tag, "PN");

	this->tag = string(tag);

}

EccBinary::~EccBinary() {
	// TODO Auto-generated destructor stub

	delete this->orb;
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


double check_LSO_condition(double F, double et, EccBinary *eb)
{
	// if this condition is negative the Last Stable Orbit (LSO) has not been breached
	double condition;
	double m = eb->get_m();

	condition = (1. + et)/(6. + 2.*et);
	condition = F - pow(condition, 1.5)/(PI2*m);

	return condition;
}







//} /* namespace eb */


