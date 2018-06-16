/*
 * EccOrbit.cpp
 *
 *  Created on: Jun 5, 2018
 *      Author: Travis Robson
 */

#include <cmath>
#include <fstream>
#include <iostream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>

#include "EccBinary.h"
#include "Constants.h"

using namespace std;

Orbit::Orbit() {
	// TODO Auto-generated constructor stub

}

Orbit::~Orbit() {
	// TODO Auto-generated destructor stub
}

double get_etdot_0PN(double et)
{
	double edot, esq;

	esq = et*et;

	edot = pow(1. - esq, -2.5)*(304./15. + 121./15.*esq);

	return  edot;
}

double get_etdot_1PN(double et, EccBinary *eb)
{
	double edot, esq;
	double eta = eb->get_eta();

	esq = et*et;

	edot  = -939./35. - 4084./45.*eta + esq*(29917./105. - 7753./30.*eta);
	edot *= pow(1. - esq, -3.5);

	return  edot;
}

double get_etdot(double et, double x, EccBinary *eb)
{
	int PN = eb->get_PN();
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
	ndot *= pow(1. - esq, -4.5);

	return ndot;
}

double get_ndot(double et, double x, EccBinary *eb)
{
	int PN = eb->get_PN();
	double ndot;

	double m   = eb->get_m();
	double eta = eb->get_eta();

	ndot = get_ndot_0PN(et);

	if (PN > 0) ndot += x*get_ndot_1PN(et, eb);

	ndot *= pow(x, 5.5)*eta/m/m*ndot;

	return ndot;
}


double get_x(double et, double n, EccBinary *eb)
{
	int PN = eb->get_PN();
	double x, zeta;
	double m   = eb->get_m();
	// double eta = eb->get_eta(); Needed at 2PN I believe

	zeta = pow(m*n, 2./3.);

	x = 1.;

	if (PN > 0) x += 2./(1. - et*et)*zeta;

	x *= zeta;

	return x;
}

int jac_EccOrbit(double t, const double y[], double *dfdy, double dfdt[], void *params)
{
	(void)(t); // avoid unused parameter warning
//	double m00, m01, m02;
//	double m10, m11, m12;
//	double m20, m21, m22;
//
//	struct EccBinary *eb = (struct EccBinary *)params;
//	double m, eta;
//
//	double e, p, phi;
//
//	e     = y[0]; // eccentricity
//	p     = y[1]; // dimensionless semi-latus rectum
//	phi   = y[2]; // orbital phase
//
//	double eSQ  = e*e;
//	double temp = pow(1. - eSQ, 1.5); // a common term
//
//	m   = eb->m;   // total mass
//	eta = eb->eta; // symmetric mass ratio
//
//  	gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 3, 3);
//  	gsl_matrix *mat = &dfdy_mat.matrix;
	  	gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 2, 2);
	  	gsl_matrix *mat = &dfdy_mat.matrix;
//
//  	// first row
//  	m00 = sqrt(1. - eSQ)*eta*(-304. + 853.*eSQ + 726.*eSQ*eSQ)/15.*pow(p, -4.);
//  	m01 = 4.*e*temp*eta*(304. + 121*eSQ)/15.*pow(p, -5.);
//  	m02 = 0.;
//
//  	// second row
//  	m10 = 8.*e*sqrt(1. - eSQ)*eta*(2. + 7.*eSQ)*pow(p, -3.);
//  	m11 = 24.*temp*eta*(8. + 7.*eSQ)/5.*pow(p, -4.);
//  	m12 = 0.;
//
//    // eighth row
//  	m20 = 2.*cos(phi)*(1. + e*cos(phi))*pow(p, -1.5);
//  	m21 = -3./2.*(1. + e*cos(phi))*(1. + e*cos(phi))*pow(p, -2.5);
//  	m22 = -2.*e*(1. + e*cos(phi))*sin(phi)*pow(p, -1.5);
//
//	// first row
//	gsl_matrix_set(mat, 0, 0, m00);
//	gsl_matrix_set(mat, 0, 1, m01);
//	gsl_matrix_set(mat, 0, 2, m02);
//
//	// second row
//	gsl_matrix_set(mat, 1, 0, m10);
//	gsl_matrix_set(mat, 1, 1, m11);
// 	gsl_matrix_set(mat, 1, 2, m12);
//
// 	// third row
//	gsl_matrix_set(mat, 2, 0, m20);
//	gsl_matrix_set(mat, 2, 1, m21);
//	gsl_matrix_set(mat, 2, 2, m22);

		// first row
		gsl_matrix_set(mat, 0, 0, 0.0);
		gsl_matrix_set(mat, 0, 1, 0.0);

		// second row
		gsl_matrix_set(mat, 1, 0, 0.0);
		gsl_matrix_set(mat, 1, 1, 0.0);

	dfdt[0] = 0.0;
	dfdt[1] = 0.0;
	dfdt[2] = 0.0;
	dfdt[3] = 0.0;

	return GSL_SUCCESS;
}

int func_EccOrbit(double t, const double y[], double f[], void *params)
{
	(void)(t); // avoid unused parameter warning
	double m, eta;
	double et, n;//, phi;
	double x;

	EccBinary *eb = (EccBinary *)params;

	m     = eb->get_m();   // total mass
	eta   = eb->get_eta(); // symmetric mass ratio


	// time derivatives of quantities (dimensionless time!, t/m)
	double etdot, ndot;//, phiDOT;

	et  = y[0]; // eccentricity
	n   = y[1]; // mean motion
//	phi = y[2]; // orbital phase

	x = get_x(et, n, eb);
	cout << "x..... " << x << endl;
	etdot = get_etdot(et, x, eb);
	ndot  = get_ndot(et, x, eb);
	cout << ndot << endl;

//	pDOT   = -64./5.*eta*pow(p, -3.)*temp*(1. + 7./8.*eSQ);
//	phiDOT = (1. + e*cos(phi))*(1. + e*cos(phi))*pow(p, -1.5);

	f[0] = etdot;
	f[1] = ndot;
//	f[2] = phiDOT;
//
	return GSL_SUCCESS;
}

void evolve_EccOrbit(EccBinary *eb)
{
	int status;

	long i;

	double t1 = 1.0e20;
	double m  = eb->get_m();
	double t, dt;
	double F, et, n;
	double LSO_condition;
//	double y[3] = {eb->e0, eb->p0, eb->phi0};
	double y[2] = {eb->get_et0(), PI2*eb->get_F0()};

	string filename;

	ofstream file;

	filename  = string("Data/orbit_soln_");
	filename += eb->get_tag();
	filename += string(".dat");

	cout << "File created......... " << filename << endl;

	 file.open (filename);

 	const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk8pd; //gsl_odeiv2_step_bsimp;

//	gsl_odeiv2_step    *s   = gsl_odeiv2_step_alloc(T, 3);
 	gsl_odeiv2_step    *s   = gsl_odeiv2_step_alloc(T, 2);
	gsl_odeiv2_control *c   = gsl_odeiv2_control_y_new(1.0e-10, 1.0e-10);
//	gsl_odeiv2_evolve  *evo = gsl_odeiv2_evolve_alloc(3);
	gsl_odeiv2_evolve  *evo = gsl_odeiv2_evolve_alloc(2);


// 	gsl_odeiv2_system sys = {func_ecc_0PN, jac_ecc_0PN, 3, eb};
	gsl_odeiv2_system sys = {func_EccOrbit, jac_EccOrbit, 2, eb};

 	// Initial quantities
 	t  = 0.;
 	et = y[0];
 	n  = y[1];
 	F  = n/PI2;
 	i  = 0;     // number of samples

	// Evolve the system
 	LSO_condition = check_LSO_condition(F, et, eb);
 	cout << "LSO condition.... " << LSO_condition << endl;
 	while (LSO_condition < 0. && i < 100)
	{
		et = y[0];
		n  = y[1];


//		fprintf(file, "%.12g %.12g %.12g %.12g\n", t, y[0], y[1], y[2]);
		file << t << " " <<  y[0] << " " << y[1] << endl;

 		dt = 0.2/F; ///eb->get_m(); // estimate of a good step size

 		status = gsl_odeiv2_evolve_apply(evo, c, s, &sys, &t, t1, &dt, y);

		if (status != GSL_SUCCESS)
		{
			//printf ("error, return value=%d\n", status);
			cout << "error, return value=" << status << endl;
			break;
		}

		i++;
		F  = n/PI2;  cout << "F.... " << F << "    i...." << i << endl;
		LSO_condition = check_LSO_condition(F, et, eb);
		cout << "LSO condition.... " << LSO_condition << endl;
	}
//	eb->FLSO = Forb;
//	eb->eLSO = e;

 	gsl_odeiv2_evolve_free(evo);
 	gsl_odeiv2_control_free(c);
 	gsl_odeiv2_step_free(s);

//	eb->orbit->N = i;

	file.close();

	return;
}




