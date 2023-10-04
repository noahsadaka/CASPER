// Standard Imports
#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <vector>
#include <string>

// SPICE toolkit
#include "SpiceUsr.h"

// User-created headers
#include "SpiceBodies.h"
#include "ephem_utils.h"

//Pybind and Boost
#include "pybind11/pybind11.h"
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std; 
namespace py = pybind11;

class nbody{

};

PYBIND11_MODULE(casper, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    py::class_<SpiceBody>(m, "SpiceBody")
	    .def(py::init<const string &, const int &, const double &>())
	    .def_readonly("name", &SpiceBody::name)
	    .def_readonly("ID", &SpiceBody::ID)
	    .def_readonly("mu", &SpiceBody::mu);

    m.def("load_kernels", &load_kernels);

    m.def("str2et", &str2et);

    py::class_<nbody>(m, "nbody")
	    .def(py::init<>());
}


