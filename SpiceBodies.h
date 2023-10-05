#ifndef CASPER_SPICE_BODIES_H
#define CASPER_SPICE_BODIES_H

#include <string>
using namespace std;

class SpiceBody{
	// Class to hold SPICE body parameters
	public:
		string name; // Human-readeable name. Not used in spice calls
		int ID; // SPICE ID. Used in spice calls
		double mu; // Gravitational parameter of body
		SpiceBody(const string n, const int I, const double m){
			name = n;
			ID = I;
			mu = m;
		};
};

#endif
