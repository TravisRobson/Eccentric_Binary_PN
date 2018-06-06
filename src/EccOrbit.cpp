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

double Orbit::get_etdot_0PN(double et)
{
	double edot, esq;

	esq = et*et;

	edot = pow(1. - esq, -2.5)*(304./15. + 121./15.*esq);

	return  edot;
}

double Orbit::get_etdot_1PN(double et, EccBinary *eb)
{
	double edot, esq;
	double eta = eb->get_eta();

	esq = et*et;

	edot  = -939./35. - 4084./45.*eta + esq*(29917./105. - 7753./30.*eta);
	edot *= pow(1. - esq, -3.5);

	return  edot;
}

double Orbit::get_etdot(double et, double x, EccBinary *eb, int PN)
{
	double edot;

	double m   = eb->get_m();
	double eta = eb->get_eta();

	edot = get_etdot_0PN(et);

	if (PN > 0) edot += x*get_etdot_1PN(et, eb);

	edot *= -et*x*x*x*x*eta/m;

	return edot;
}



double get_ndot_0PN(double et)
{
	double ndot;
	double esq = et*et;

	ndot  = 96./5. + 292./5.*esq + 37./5.*esq*esq;
	ndot *= pow(1. - esq, -3.5);

	return ndot;
}

double get_ndot_1PN(double et, EccBinary *eb)
{
	double ndot;
	double eta = eb->get_eta();
	double esq = et*et;

	ndot  = -4846./35. - 264./5.*eta + esq*(5001./35. - 570.*eta);
	ndot += esq*esq*(2489./4. - 5061./10.*eta) + esq*esq*esq*(11717./280. - 148./5.*eta);
	ndot *= pow(1. -esq, -4.5);

	return ndot;
}

double get_ndot(double et, double x, EccBinary *eb, int PN)
{
	double ndot;

	double m   = eb->get_m();
	double eta = eb->get_eta();

	ndot = get_ndot_0PN(et);

	if (PN > 0) ndot += x*get_ndot_1PN(et, eb);

	ndot *= pow(x, 5.5)*eta/m/m*ndot;

	return ndot;
}




