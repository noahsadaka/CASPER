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

class cr3bp_system{
	public:
		double m_star;
		double l_star;
		double t_star;
	cr3bp_system(const double m, const double l, const double t){
		m_star = m;
		l_star = l;
		t_star = t;
	};
};

class nbody{

};

PYBIND11_MODULE(casper, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    py::class_<SpiceBody>(m, "SpiceBody")
	    .def(py::init<const string &, const int &, const double &>())
	    .def_readonly("name", &SpiceBody::name)
	    .def_readonly("ID", &SpiceBody::ID)
	    .def_readonly("mu", &SpiceBody::mu);

    //m.def("load_kernels", &load_kernels);
    m.def("str2et", &str2et);

    py::class_<nbody>(m, "nbody")
	    .def(py::init<>());

    py::class_<cr3bp_system>(m, "cr3bp_system")
	    .def(py::init<const double &, const double &, const double &>())
	    .def_readonly("m_star", &cr3bp_system::m_star)
	    .def_readonly("l_star", &cr3bp_system::l_star)
	    .def_readonly("t_star", &cr3bp_system::t_star);
}


