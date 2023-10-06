#ifndef CASPER_SPICE_BODIES_H
#define CASPER_SPICE_BODIES_H

#include <string>

#include "ephem_utils.h"

using namespace std;

class SpiceEpoch{
    public:
        // Constructors, can handle input as ephemeris time float
        // or as UTC time string with whatever format str2et_c accepts
        SpiceEpoch(SpiceDouble et){
            load_kernels();
            et_epoch = et;
            utc_epoch = et2utc(et);
        };
        SpiceEpoch(string utc){
            load_kernels();
            et_epoch = str2et(utc);
            utc_epoch = utc;
        };

        //Attributes
        string utc_epoch;
        SpiceDouble et_epoch;
};

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
