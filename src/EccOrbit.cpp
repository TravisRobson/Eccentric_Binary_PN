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

	edot = pow(1. - esq, -2.5)*(304. + 121.*esq)/15.;

	return  edot;
}

double get_etdot_1PN(double et, EccBinary *eb)
{
	double edot, esq;
	double eta = eb->get_eta();

	esq = et*et;

	edot  = 340968. - 228704.*eta + esq*(880632. - 651252.*eta);
	edot += esq*esq*(125361. - 93184*eta);
	edot /= 2520.;
	edot *= pow(1. - esq, -3.5);

	return  edot;
}


double get_etdot(double et, double zeta, EccBinary *eb)
{
	int PN = eb->get_PN();
	double edot;

	double eta = eb->get_eta();

	edot = get_etdot_0PN(et);

	if (PN > 0) edot += pow(zeta, 2./3.)*get_etdot_1PN(et, eb);

	edot *= -et*pow(zeta, 8./3.)*eta; // 1./m neglected to make it dimensionless time

	return edot;
}



double get_ndot_0PN(double et)
{
	double ndot;
	double esq = et*et;

	ndot  = (96. + 292.*esq + 37.*esq*esq)/5.;
	ndot *= pow(1. - esq, -3.5);

	return ndot;
}

double get_ndot_1PN(double et, EccBinary *eb)
{
	double ndot;
	double eta = eb->get_eta();
	double esq = et*et;

	ndot  = 20368. - 14784*eta + esq*(219880. - 159600.*eta);
	ndot += esq*esq*(197022 - 141708.*eta) + esq*esq*esq*(11717. - 8288.*eta);
	ndot /= 280.;
	ndot *= pow(1. - esq, -4.5);

	return ndot;
}

double get_ndot(double et, double zeta, EccBinary *eb)
{
	int PN = eb->get_PN();
	double ndot;

	double eta = eb->get_eta();

	ndot = get_ndot_0PN(et);

	if (PN > 0) ndot += pow(zeta, 2./3.)*get_ndot_1PN(et, eb);

	ndot *= pow(zeta, 11./3.)*eta; // 1./m**2 neglected to make dimensionless (mn) and dimensionless time

	return ndot;
}


double get_x(double et, double n, EccBinary *eb)
{
	int PN = eb->get_PN();
	double x, zeta;
	double m   = eb->get_m();

	zeta = pow(m*n, 2./3.);

	x = 1.;

	if (PN > 0) x += 2./(1. - et*et)*zeta;

	x *= zeta;

	return x;
}

double get_zeta(double n, EccBinary *eb)
{
	return eb->get_m()*n;
}

double get_k(double zeta, double et)
{
	double k;

	k = 3.*pow(zeta, 2./3.)/(1. - et*et);

	return k;
}

double get_phidot(double et, double n, double phi, EccBinary *eb)
{
	int PN = eb->get_PN();
	double phidot;
	double k, zeta, nu;
	double term;
	double eta = eb->get_eta();
	double m = eb->get_m();

	// get zeta
	zeta = get_zeta(n, eb);

	// get k
	if (PN > 0) k = get_k(zeta, et);
	else k = 1.;

	// get true anomaly
	nu = phi/(1. + k);

	// get 0 pN term
	phidot  = n*(1. + et*cos(nu))*(1. + et*cos(nu))/pow(1. - et*et, 1.5);
	phidot += n;
	if (PN > 0) phidot += k*n;


	// get 1 pN term
	if (PN > 0)
	{
		term  = -(4. + 5.*et*et)*cos(nu) + et*(-8. + et*cos(nu));
		term *= (4. - eta)*et/pow(1. -et*et, 2.5);

		phidot += pow(zeta, 2./3.)*term*n;
	}


	return phidot*m; // multiply by m for dimensionless time
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
	gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 3, 3);
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
	gsl_matrix_set(mat, 0, 0, 1.0);
	gsl_matrix_set(mat, 0, 1, 0.0);
	gsl_matrix_set(mat, 0, 2, 0.0);

	// second row
	gsl_matrix_set(mat, 1, 0, 0.0);
	gsl_matrix_set(mat, 1, 1, 1.0);
	gsl_matrix_set(mat, 1, 2, 0.0);

	// third row
	gsl_matrix_set(mat, 2, 0, 0.0);
	gsl_matrix_set(mat, 2, 1, 0.0);
	gsl_matrix_set(mat, 2, 2, 1.0);


	dfdt[0] = 0.0;
	dfdt[1] = 0.0;
	dfdt[2] = 0.0;
//	dfdt[3] = 0.0;

	return GSL_SUCCESS;
}

int func_EccOrbit(double t, const double y[], double f[], void *params)
{
	(void)(t); // avoid unused parameter warning
	double m, eta;
	double et, n, phi;
	double zeta;

	EccBinary *eb = (EccBinary *)params;

	m     = eb->get_m();   // total mass
	eta   = eb->get_eta(); // symmetric mass ratio

	// time derivatives of quantities (dimensionless time!, t/m)
	double etdot, ndot, phidot;

	et  = y[0];   // eccentricity
	n   = y[1]/m; // mean motion
	phi = y[2];   // orbital phase

	zeta = get_zeta(n, eb);

	etdot = get_etdot(et, zeta, eb);
	ndot  = get_ndot(et, zeta, eb);

	phidot = get_phidot(et, n, phi, eb);

	f[0] = etdot;
	f[1] = ndot;
	f[2] = phidot;

	return GSL_SUCCESS;
}

void evolve_EccOrbit(EccBinary *eb)
{
	int status;

	long i;

	double t1 = 1.0e20;
	double m  = eb->get_m();
	double t, dt;
	double F, et, n, phi;
	double LSO_condition;
	double y[3] = {eb->get_et0(), PI2*eb->get_F0()*m, eb->get_phi0()};
//	double y[2] = {eb->get_et0(), PI2*eb->get_F0()*m};

	string filename;

	ofstream file;

	filename  = string("Data/orbit_soln_");
	filename += eb->get_tag();
	filename += string(".dat");

	cout << "File created......... " << filename << endl;

	file.open (filename);

 	const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk8pd; //gsl_odeiv2_step_bsimp;
 	//const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rkf45;

	gsl_odeiv2_step    *s   = gsl_odeiv2_step_alloc(T, 3);
// 	gsl_odeiv2_step    *s   = gsl_odeiv2_step_alloc(T, 2);
	gsl_odeiv2_control *c   = gsl_odeiv2_control_y_new(1.0e-16, 1.0e-16);
	gsl_odeiv2_evolve  *evo = gsl_odeiv2_evolve_alloc(3);
//	gsl_odeiv2_evolve  *evo = gsl_odeiv2_evolve_alloc(2);


 	gsl_odeiv2_system sys = {func_EccOrbit, jac_EccOrbit, 3, eb};
//	gsl_odeiv2_system sys = {func_EccOrbit, jac_EccOrbit, 2, eb};

 	// Initial quantities
 	t   = 0.;
 	et  = y[0];
 	n   = y[1]/m;
 	phi = y[2];
 	F   = n/PI2;
 	i   = 0;     // number of samples

	// Evolve the system
 	LSO_condition = check_LSO_condition(F, et, eb);

 	while (LSO_condition < 0.)// && i<214)
	{
		et  = y[0];
		n   = y[1]/m;
		F   = n/PI2;
		phi = y[2];

		file << t*m << " " <<  et << " " << F << " " << phi << endl;

 		dt = 0.1/F/m; // estimate of a good step size

 		status = gsl_odeiv2_evolve_apply(evo, c, s, &sys, &t, t1, &dt, y);

		if (status != GSL_SUCCESS)
		{
			cout << "error, return value=" << status << endl;
			break;
		}

		i++;
		LSO_condition = check_LSO_condition(F, et, eb);
	}
	eb->FLSO = F;
	eb->etLSO = et;

 	gsl_odeiv2_evolve_free(evo);
 	gsl_odeiv2_control_free(c);
 	gsl_odeiv2_step_free(s);

	//eb->orbit->N = i;

	file.close();

	return;
}




