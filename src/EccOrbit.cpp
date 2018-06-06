/*
 * EccOrbit.cpp
 *
 *  Created on: Jun 5, 2018
 *      Author: travisrobson
 */

#include <cmath>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>

#include "EccBinary.h"
#include "Constants.h"

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
	ndot *= pow(1. -esq, -4.5);

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

	zeta = pow(m*n, 2./3.); ;

	x = 1.;

	if (PN > 0) x += 2./(1. - et*et)*zeta;

	x *= zeta;

	return x;
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

	etdot = get_etdot(et, x, eb);
	ndot  = get_ndot(et, x, eb);

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
	double t, dt, Forb, e, p, LSO_condition;
//	double y[3] = {eb->e0, eb->p0, eb->phi0};
	double y[2] = {eb->get_et0(), PI2*eb->get_F0()};
//
//	char buf[512];
//
//	FILE *file;
//
//	strcpy(buf, "orbit_soln_");
//	strcat(buf, eb->tag);
//	strcat(buf, ".dat");
//
//	file = fopen(buf, "w");
//
// 	const gsl_odeiv2_step_type * T = gsl_odeiv2_step_bsimp;
//
//	gsl_odeiv2_step    *s   = gsl_odeiv2_step_alloc(T, 3);
//	gsl_odeiv2_control *c   = gsl_odeiv2_control_y_new(1.0e-16, 1.0e-16);
//	gsl_odeiv2_evolve  *evo = gsl_odeiv2_evolve_alloc(3);
//
// 	gsl_odeiv2_system sys = {func_ecc_0PN, jac_ecc_0PN, 3, eb};
//
// 	/* ---- Initial quantities ---- */
// 	t = 0.; // initial time
// 	e = y[0];
// 	p = y[1];
// 	Forb = pow((1. - e*e)/p, 3./2.)/(PI2*m);
// 	i = 0; // number of samples
//
// 	/* ---- Evolve the system ----*/
// 	LSO_condition = Forb - pow((1. + e)/(6. + 2.*e), 1.5)/(PI2*m);
// 	while (LSO_condition < 0.)
//	{
//		e    = y[0];
//		p    = y[1];
//		Forb = pow((1. - e*e)/p, 3./2.)/(PI2*m);
//
//		fprintf(file, "%.12g %.12g %.12g %.12g\n", t, y[0], y[1], y[2]);
//
// 		dt = 0.2/Forb/m; // estimate of a good step size
//
// 		status = gsl_odeiv2_evolve_apply(evo, c, s, &sys, &t, t1, &dt, y);
//
//		if (status != GSL_SUCCESS)
//		{
//			printf ("error, return value=%d\n", status);
//			break;
//		}
//
//		i++;
//		LSO_condition = Forb - pow((1. + e)/(6. + 2.*e), 1.5)/(PI2*m);
//	}
//	eb->FLSO = Forb;
//	eb->eLSO = e;
//
// 	gsl_odeiv2_evolve_free(evo);
// 	gsl_odeiv2_control_free(c);
// 	gsl_odeiv2_step_free(s);
//
//	eb->orbit->N = i;
//
//	fclose(file);

	return;
}




