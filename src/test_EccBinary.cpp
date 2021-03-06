/*
 * test_EccBinary.c
 *
 *  Created on: Jun 5, 2018
 *      Author: Travis Robson
 */

#include <iostream>
#include <time.h>

//#include <string>
//#include <sstream>

//#include <cstdlib>

using namespace std;

#include "EccBinary.h"
#include "Constants.h"

int main(int argc, char **argv)
{
	int PN;

	double m1,  m2;
	double et0, F0, phi0;

	cout << "\n==================================================\n" << endl;

	PN = 1;
	m1  = 10.0*TSUN;
	m2  = 10.0*TSUN;
	et0 = 0.5;
	F0  = 5.0;
	phi0 = 0.0;

	EccBinary *eb = new EccBinary(m1, m2, et0, F0, phi0, PN);

	cout << "m............ " << eb->get_m()/TSUN << " MSUN"  << endl;
	cout << "eta.......... " << eb->get_eta()         	   << endl;
	cout << "Mc........... " << eb->get_Mc()/TSUN << " MSUN" << endl << endl;

	cout << "N for soln... " << eb->orb->get_N() << endl;

	cout << "tag.......... " << eb->get_tag() << endl;

	evolve_EccOrbit(eb);

	cout << "\n==================================================\n" << endl;

	delete eb;

	return 0;
}
