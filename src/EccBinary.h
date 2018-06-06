/*
 * EccBinary.h
 *
 *  Created on: Jun 5, 2018
 *      Author: Travis Robson
 */

#ifndef ECCBINARY_H_
#define ECCBINARY_H_

#include <cstdlib>

#include <gsl/gsl_spline.h>

class Orbit;
class EccBinary;

class Orbit {
private:
	long N;				 // Number of samples for the solution (NOT number of detector samples!!!)

	double *et, *n, *phi; // eccentricity, mean motion, orbital phase pointers
	double *t;			  // time associated with the above samples

	// interpolants and accelerators for the orbital parameters
	gsl_interp_accel *e_acc;
	gsl_spline *e_spline;

	gsl_interp_accel *n_acc;
	gsl_spline *n_spline;

	gsl_interp_accel *phi_acc;
	gsl_spline *phi_spline;

public:
	Orbit();
	virtual ~Orbit();

	long get_N() { return N; }
	void set_N(long no) { N = no; } // This is dummy TODO

	double get_etdot(double et, double x, EccBinary *eb, int PN);
	double get_etdot_0PN(double et);
	double get_etdot_1PN(double et, EccBinary *eb);

	double get_ndot(double et, double x, EccBinary *eb, int PN);
	double get_ndot_0PN(double et);
	double get_ndot_1PN(double et, EccBinary *eb);

};

//namespace eb {

class EccBinary {

private:
	double m1, m2;	       // component masses
	double m, mu, eta, Mc; // total mass, reduced mass, symmetric mass ratio

	double et0, F0; // initial eccentricity, initial mean orbital frequency (related to mean motion)



public:
	Orbit *orb;

	EccBinary(double m1, double m2, double et0, double F0);
	virtual ~EccBinary();

	// Member function declarations
	double get_m1() {
		return m1;
	}
	double get_m2() {
		return m2;
	}
	double get_m() {
		return m;
	}
	double get_mu() {
		return mu;
	}
	double get_eta() {
		return eta;
	}
	double get_Mc() {
		return Mc;
	}

	void calc_m();
	void calc_mu();
	void calc_eta();
	void calc_Mc();
	void calc_mass_parameters();
};

//} /* namespace eb */



#endif /* ECCBINARY_H_ */

