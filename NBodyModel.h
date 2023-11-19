#ifndef CASPER_NBODYMODEL_H
#define CASPER_NBODYMODEL_H
// Standard Imports
#include <iostream>

// Spice Bodies
#include "SpiceClasses.h"
 
using namespace std;


class NBodyModel{
    public:
        cr3bp_system sys; // CR3BP system for nondimensionalization
        SpiceBody central_body; // Central body of J2000 frame
        std::vector<SpiceBody> perturbing_bodies; // List of perturbing bodies

        NBodyModel(const double m, const double l, SpiceBody central, std::vector<SpiceBody> perturbing);
};

NBodyModel::NBodyModel(const double m, const double l, SpiceBody central, std::vector<SpiceBody> perturbing):
    sys(m, l),
    central_body(central),
    perturbing_bodies(perturbing)
{
};

#endif
