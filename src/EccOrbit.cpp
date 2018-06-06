/*
 * EccOrbit.cpp
 *
 *  Created on: Jun 5, 2018
 *      Author: travisrobson
 */

#include <cmath>

#include "EccBinary.h"

Orbit::Orbit() {
	// TODO Auto-generated constructor stub

}

Orbit::~Orbit() {
	// TODO Auto-generated destructor stub
}

double Orbit::get_edot_0PN(double et)
{
	double edot, esq;

	esq = et*et;

	edot = pow(1. - esq, -2.5)*(304./15. + 121./15.*esq);

	return  edot;
}

double Orbit::get_edot_1PN(double et, double eta)
{
	double edot, esq;

	esq = et*et;

	edot  = -939./35. - 4084./45.*eta + esq*(29917./105. - 7753./30.*eta);
	edot *= pow(1. - esq, -3.5);

	return  edot;
}
