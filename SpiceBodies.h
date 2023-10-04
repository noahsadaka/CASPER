#ifndef CASPER_SPICE_BODIES_H
#define CASPER_SPICE_BODIES_H

#include <string>
using namespace std;

class SpiceBody{
	public:
		string name;
		int ID;
		double mu;
		SpiceBody(const string n, const int I, const double m){
			name = n;
			ID = I;
			mu = m;
		};
};

#endif
