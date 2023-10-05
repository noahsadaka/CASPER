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
#ifdef PYTHON_COMPILE
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#endif
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std; 
using namespace boost::numeric::odeint;

#ifdef PYTHON_COMPILE
namespace py = pybind11;
#endif

typedef std::vector<double> state_type;

class PropObserver{
	public:
		std::vector<std::vector<double>> x; // integration states
		state_type t; // integration times

		void operator()(const state_type &x_curr, const double t_curr){
			t.push_back(t_curr);
			x.push_back(x_curr);
		};
};


class cr3bp_system{
	public:
		double m_star;
		double l_star;
		double t_star;
	cr3bp_system(const double m, const double l){
		const double G ((6.67430e-11/pow(1000,3)));
		m_star = m;
		l_star = l;
		t_star = 1/pow(m*G/pow(l,3), 0.5);
	};
};

class NBODY{
	public:
		
		// Attributes
		state_type IC; // initial condition, 42-vector, non-dim
		SpiceDouble base_epoch_nd; // base epoch, in non-dim seconds ephemeris time
		cr3bp_system sys; // CR3BP object used for non-dimensionalization
		SpiceBody central_body; // Body at center of inertial J2000 frame
		std::vector<SpiceBody> perturbing_bodies; // Perturbing bodies
		
		// Constructor Declaration
		NBODY(state_type &ic, string str_epoch, const double m, const double l,
			SpiceBody central, std::vector<SpiceBody> perturbing);

		// EOM function declaration
		void EOM_STM(state_type &x, state_type &dx, const double t);

		// Primary acceleration computation function declaration
		void get_primary_acceleration(const state_type &x, state_type &acc, state_type &A_subset);

		// Perturbing acceleration computation function declaration
		void get_perturbing_acceleration(const double t, SpiceBody pert_body, int central_body_ID, cr3bp_system sys, const state_type &state, state_type &acc, state_type &A_subset);

		// A matrix creation function declaration
		std::vector<double> build_A_matrix(state_type &A_sub);

		// propagator function
		void propagate(double TOF, double step_size, double rt, double at, PropObserver o);
};

// Define NBODY class Constructor
NBODY::NBODY(state_type &ic, string str_epoch, const double m, const double l,
	SpiceBody central, std::vector<SpiceBody> perturbing):
	// member definition, necessary for the cr3bp class
	IC(ic),
	sys(m, l),
	central_body(central),
	perturbing_bodies(perturbing)
	{
	// actual function here
	SpiceDouble base_epoch;
	base_epoch = str2et(str_epoch);
	// non-dimensionalize the epoch
	base_epoch_nd = base_epoch/sys.t_star;

	// make sure kernels are loaded
	load_kernels(); 
};

// EOMs with STM propagation
void NBODY::EOM_STM(state_type &state, state_type &d_state, const double t){
	state_type acc(3,0);
  	state_type A_subset (9,0);
	
	// get acceleration and A matrix contribution from central body
	get_primary_acceleration(state, acc, A_subset);	

	// iterate through perturbing body list to get acceleration and A matrix contributions
	for (int i = 0; i < perturbing_bodies.size(); i++)
	{
		get_perturbing_acceleration(t, perturbing_bodies[i], central_body.ID, sys, state, acc, A_subset);
	}

	// Create the 6x6 A matrix from the 3x3 subset 
	state_type A;
	A = build_A_matrix(A_subset);

	// Velocity terms
	d_state[0] = state[3];
	d_state[1] = state[4];
	d_state[2] = state[5];
	d_state[3] = acc[0];
	d_state[4] = acc[1];
	d_state[5] = acc[2];
	
	// Do the phi_dot = A * phi computation. Took this from Nick and RJ, thanks Nick and RJ!
	int state_size = 6;
	for (int i = 0; i < state_size; i++) {
              for (int j = 0; j < state_size; j++) {
                  unsigned int current_index = state_size + i * state_size + j;

                  d_state[current_index] = 0.0;

                  for (int k = 0; k < state_size; k++) {
                      d_state[current_index] += A[i * state_size + k] * state[state_size + state_size * k + j];
                  }
              }
          }
     };

// Primary Acceleration
void NBODY::get_primary_acceleration(const state_type &state, state_type &acc, state_type &A_subset){
	const double G ((6.67430e-11/pow(1000,3)));
	const double mass (central_body.mu/G/sys.m_star);
	
	// extract only position components of state
	state_type r_sc_b {state[0], state[1], state[2]};

	const double r_sc_b_n = sqrt(r_sc_b[0]*r_sc_b[0] + r_sc_b[1]*r_sc_b[1] + r_sc_b[2]*r_sc_b[2]);
	const double r_sc_b_n_5 = pow(r_sc_b_n, 5);
	const double r_sc_b_n_3 = pow(r_sc_b_n, 3);
	
	// makes it easier to write and read
	double x, y, z;
	x = r_sc_b[0];
	y = r_sc_b[1];
	z = r_sc_b[2];
	
	// acceleration terms
	acc[0] -= mass * x / r_sc_b_n_3;
	acc[1] -= mass * y / r_sc_b_n_3;
	acc[2] -= mass * z / r_sc_b_n_3;
	
	// STM terms
	A_subset[0] += mass * (3 * pow(x, 2) / r_sc_b_n_5 - 1 / r_sc_b_n_3);
	A_subset[1] += mass * 3 * x * y / r_sc_b_n_5;
	A_subset[2] += mass * 3 * x * z / r_sc_b_n_5;
	A_subset[3] += mass * 3 * x * y / r_sc_b_n_5;
	A_subset[4] += mass * (3 * pow(y, 2) / r_sc_b_n_5 - 1 / r_sc_b_n_3);
	A_subset[5] += mass * 3 * y * z / r_sc_b_n_5;
	A_subset[6] += mass * 3 * x * z / r_sc_b_n_5;
	A_subset[7] += mass * 3 * y * z / r_sc_b_n_5;
	A_subset[8] += mass * (3 * pow(z, 2) / r_sc_b_n_5 - 1 / r_sc_b_n_3);
};

// Perturbing acceleration
void NBODY::get_perturbing_acceleration(const double t, SpiceBody pert_body,
	int central_body_ID, cr3bp_system sys, const state_type &state,
	state_type &acc, state_type &A_subset){
	
	// variables
	const double G ((6.67430e-11/pow(1000,3)));
	const double mass (pert_body.mu/G/sys.m_star);
	const double t_star (sys.t_star);
	SpiceDouble body_position [3], lighttimes;
	SpiceDouble epoch_dim (t*t_star);
	
	// get position of perturbing body WRT central body	
	spkpos_c(to_string(pert_body.ID).c_str(), epoch_dim, "J2000", "NONE",
		 to_string(central_body_ID).c_str(), body_position, &lighttimes);
	//cout<< body_position[0]<<endl;
	// non dimensionalize planet position
	body_position[0] = body_position[0]/sys.l_star;
	body_position[1] = body_position[1]/sys.l_star;
	body_position[2] = body_position[2]/sys.l_star;
	
	state_type r_sc_body, r_body_obs;
	double r_sc_b_n, r_sc_b_n_3, r_sc_b_n_5, x, y, z;
	double x_o, y_o, z_o;
	double r_b_o_n, r_b_o_n_3;
	
	// Vector pointing from perturbing body to spacecraft
	r_sc_body = {state[0] - body_position[0], state[1] - body_position[1], state[2] - body_position[2]};
	// Vector pointing from central body to perturbing body
	r_body_obs = {-1 * body_position[0],-1 * body_position[1], -1 * body_position[2]};
	// norm of spacecraft to central body vector
	r_sc_b_n = sqrt(pow(r_sc_body[0],2) + pow(r_sc_body[1],2)+ pow(r_sc_body[2],2));
	r_sc_b_n_3 = pow(r_sc_b_n,3);
	r_sc_b_n_5 = pow(r_sc_b_n,5);
	// norm of perturbing to central body vector
	r_b_o_n = sqrt(pow(r_body_obs[0],2) + pow(r_body_obs[1],2)+ pow(r_body_obs[2],2));
	r_b_o_n_3 = pow(r_b_o_n,3);

	// states of r_sc_body vector
	x = r_sc_body[0];
	y = r_sc_body[1];
	z = r_sc_body[2];
	// states of r_body_obs vector
	x_o = r_body_obs[0];
	y_o = r_body_obs[1];
	z_o = r_body_obs[2];
	// add perturbing accelerations to acceleration vector
	acc[0] -= mass * (x / r_sc_b_n_3 - x_o / r_b_o_n_3);
	acc[1] -= mass * (y / r_sc_b_n_3 - y_o / r_b_o_n_3);
	acc[2] -= mass * (z / r_sc_b_n_3 - z_o / r_b_o_n_3);
	// Add perturbing contributions to A matrix
	A_subset[0] += mass * (3 * pow(x, 2) / r_sc_b_n_5 - 1 / r_sc_b_n_3);
	A_subset[1] += mass * 3 * x * y / r_sc_b_n_5;
	A_subset[2] += mass * 3 * x * z /r_sc_b_n_5;
	A_subset[3] += mass * 3 * x * y / r_sc_b_n_5;
	A_subset[4] += mass * (3 * pow(y, 2) / r_sc_b_n_5 - 1 / r_sc_b_n_3);
	A_subset[5] += mass * (3 * y * z / r_sc_b_n_5);
	A_subset[6] += mass * 3 * x * z /r_sc_b_n_5;
	A_subset[7] += mass * (3 * y * z / r_sc_b_n_5);
	A_subset[8] += mass * (3 * pow(z, 2) / r_sc_b_n_5 - 1 / r_sc_b_n_3);
};

std::vector<double> NBODY::build_A_matrix(state_type &A_sub){
	std::vector<double> A(36,0);
	// Zero matrix block
	A[3] = 1;
	A[10] = 1;
	A[17] = 1;
	// Terms from bodies block
	A[18] = A_sub[0];
	A[19] = A_sub[1];
	A[20] = A_sub[2];
	A[24] = A_sub[3];
	A[25] = A_sub[4];
	A[26] = A_sub[5];
	A[30] = A_sub[6];
	A[31] = A_sub[7];
	A[32] = A_sub[8];

	return A;
};

void NBODY::propagate(double t_end, double step_size, double rtol, double atol, PropObserver o){
	namespace pl = std::placeholders;
	state_type states_and_times;
	// check propagation direction
	if (t_end < base_epoch_nd) {step_size = - step_size;}

	typedef runge_kutta_fehlberg78<state_type> rk78;

	auto stepper = make_controlled<rk78>(atol, rtol);

	size_t steps;

	steps = integrate_adaptive(stepper, std::bind(&NBODY::EOM_STM, *this, pl::_1, pl::_2, pl::_3), IC, base_epoch_nd, t_end, step_size, std::ref(o));
};

#ifdef PYTHON_COMPILE
PYBIND11_MODULE(casper, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    py::class_<SpiceBody>(m, "SpiceBody")
	    .def(py::init<const string &, const int &, const double &>())
	    .def_readonly("name", &SpiceBody::name)
	    .def_readonly("ID", &SpiceBody::ID)
	    .def_readonly("mu", &SpiceBody::mu);

    m.def("str2et", &str2et);

    py::class_<NBODY>(m, "nbody")
	    .def(py::init<state_type &, string, const double, const double, SpiceBody, std::vector<SpiceBody>>())
	    .def_readonly("base_epoch", &NBODY::base_epoch_nd)
	    .def_readonly("IC", &NBODY::IC)
	    .def("propagate", &NBODY::propagate);

    py::class_<PropObserver>(m, "PropObserver")
	    .def(py::init<>())
	    .def_readonly("t", &PropObserver::t)
	    .def_readonly("x", &PropObserver::x);

    //py::class_<cr3bp_system>(m, "cr3bp_system")
//	    .def(py::init<const double, const double>())
//	    .def_readonly("m_star", &cr3bp_system::m_star)
//	    .def_readonly("l_star", &cr3bp_system::l_star)
//	    .def_readonly("t_star", &cr3bp_system::t_star);
}
#endif

int main(){
	cout<<"testing file starts here"<<endl;

	const double G ((6.67430e-11/pow(1000,3))); //  # km3/kg/s2

	SpiceBody Earth("EARTH", 399, 3.9860043543609593e+05);
	SpiceBody Moon("MOON", 301, 4.9028000661637961e+03);
	SpiceBody Sun("SUN", 10, 1.3271244004193930e+11);
	SpiceBody Jupiter("JUPITER BARYCENTER", 5, 126712767.8578);
	//cr3bp_system crtbp_sys((Earth.mu+Moon.mu)/G, 3.8474799197904585e+05);
	//cout << setprecision(17) << crtbp_sys.t_star << endl;
	state_type IC_vector {1.05903, -0.067492, -0.103524, -0.170109, 0.0960234, -0.135279, 1,
          1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0};
	const double m_star ((Earth.mu+Moon.mu)/G);
	const double l_star (3.8474799197904585e+05);
	NBODY nbody_system(IC_vector, "May 2, 2022", m_star, l_star, Earth, {Moon});
	nbody_system.EOM_STM(nbody_system.IC, nbody_system.IC, nbody_system.base_epoch_nd);	
	PropObserver obs{};
	nbody_system.propagate(nbody_system.base_epoch_nd + 3, 1e-5, 1e-12, 1e-12, obs);
};
