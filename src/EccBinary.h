/*
 * EccBinary.h
 *
 *  Created on: Jun 5, 2018
 *      Author: Travis Robson
 */

#ifndef ECCBINARY_H_
#define ECCBINARY_H_

#include <cstdlib>

//namespace eb {

class EccBinary
{

private:
	double m1, m2;	   // component masses
	//double m, mu, eta; // total mass, reduced mass, symmetric mass ratio

	double e0, F0;     // initial eccentricity, initial mean orbital frequency

	public:
		EccBinary(double m1, double m2, double e0, double F0);
		virtual ~EccBinary();

		// Member function declarations
		double get_m1() { return m1; }
		double get_m2() { return m2; }
		//void calc_mu():

};

//} /* namespace eb */





#endif /* ECCBINARY_H_ */
