/*
 * test_EccBinary.c
 *
 *  Created on: Jun 5, 2018
 *      Author: travisrobson
 */

//#include <string>
//#include <sstream>
#include <iostream>
//#include <cstdlib>
#include <time.h>
using namespace std;

#include "EccBinary.h"

int main(int argc, char **argv)
{
	double m1, m2;
	double e0, F0;

	cout << "==================================================\n" << endl;

	m1 = 10.0*4.93e-6;
	m2 = 10.0*4.93e-6;
	e0 = 0.1;
	F0 = 5.0;

	EccBinary *eb = new EccBinary(m1, m2, e0, F0);

	cout << "m1........." << eb->get_m1()/4.93e-6 << " MSUN"<< endl;


	cout << "\n==================================================\n" << endl;

	return 0;
}
