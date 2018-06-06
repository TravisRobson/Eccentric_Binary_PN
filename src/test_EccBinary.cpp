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
	double m1, m2;
	double e0, F0;

	cout << "\n==================================================\n" << endl;

	m1 = 10.0*TSUN;
	m2 = 10.0*TSUN;
	e0 = 0.1;
	F0 = 5.0;

	EccBinary *eb = new EccBinary(m1, m2, e0, F0);

	cout << "m.......... " << eb->get_m()/TSUN << " MSUN"  << endl;
	cout << "eta........ " << eb->get_eta()         	   << endl;
	cout << "Mc......... " << eb->get_Mc()/TSUN << " MSUN" << endl;


	cout << "\n==================================================\n" << endl;

	return 0;
}
