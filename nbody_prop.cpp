// Standard Imports
#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <vector>
#include <string>

// SPICE toolkit
#include "SpiceUsr.h"

// User-created headers
#include "SpiceClasses.h"
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
		const unsigned int n_IC; // number of states, should be 6 or 42
		bool prop_STM;
		SpiceDouble base_epoch_nd; // base epoch, in non-dim seconds ephemeris time
		cr3bp_system sys; // CR3BP object used for non-dimensionalization
		SpiceBody central_body; // Body at center of inertial J2000 frame
		std::vector<SpiceBody> perturbing_bodies; // Perturbing bodies
        const double G = (6.67430e-11/pow(1000,3));
		
		// Constructor Declaration
		NBODY(state_type &ic, SpiceEpoch epoch, const double m, const double l,
			SpiceBody central, std::vector<SpiceBody> perturbing);

		// EOM function declaration
		void EOM_STM(state_type &x, state_type &dx, const double t);

		// Primary acceleration computation function declaration
		void get_primary_acceleration(const state_type &x, state_type &acc, state_type &A_subset);

		// Perturbing acceleration computation function declaration
		void get_perturbing_acceleration(const double t, SpiceBody pert_body,
				int central_body_ID, cr3bp_system sys, const state_type &state,
				state_type &acc, state_type &A_subset, state_type &sum_term);

        // partial of states wrt epoch summation term computation (eq 4.11 in ATD mathspec)
        void states_wrt_epoch_summation(state_type r_s_i, state_type r_c_i,
                const double body_mass, state_type &sum_term);

		// A matrix creation function declaration
		std::vector<double> build_A_matrix(state_type &A_sub);

		// propagator function
		void propagate(double TOF, double step_size, double rt, double at, PropObserver &o);
};

// Define NBODY class Constructor
NBODY::NBODY(state_type &ic, SpiceEpoch epoch, const double m, const double l,
	SpiceBody central, std::vector<SpiceBody> perturbing):
	// member definition, necessary for the cr3bp class
	IC(ic),
	sys(m, l),
	n_IC(IC.size()),
	central_body(central),
	perturbing_bodies(perturbing)
	{// actual function here
	// Check that len(IC) = 6 or 48
	if (n_IC != 6 && n_IC != 48){
		throw std::length_error("IC is not length 6 or 48");
	} else if (n_IC == 6){prop_STM = false;}
	else {prop_STM = true;}
    
    // non-dimensionalize the epoch
	base_epoch_nd = epoch.et_epoch/sys.t_star;

	// make sure kernels are loaded
	load_kernels(); 
};

// EOMs with STM propagation
void NBODY::EOM_STM(state_type &state, state_type &d_state, const double t){
	state_type acc(3, 0);
  	state_type A_subset (9, 0);
    state_type epoch_partial_sum_term (6, 0);
	
	// get acceleration and A matrix contribution from central body
	get_primary_acceleration(state, acc, A_subset);	

	// iterate through perturbing body list to get acceleration and A matrix contributions
	for (int i = 0; i < perturbing_bodies.size(); i++)
	{
		get_perturbing_acceleration(t, perturbing_bodies[i], central_body.ID, sys, state,
				acc, A_subset, epoch_partial_sum_term);
	}

    // Velocity terms
	d_state[0] = state[3];
	d_state[1] = state[4];
	d_state[2] = state[5];
	d_state[3] = acc[0];
	d_state[4] = acc[1];
	d_state[5] = acc[2];

	if (prop_STM) {
		// Create the 6x6 A matrix from the 3x3 subset 
		state_type A;
		A = build_A_matrix(A_subset);
		// Do the phi_dot = A * phi computation. Took this from Nick and RJ, thanks Nick and RJ!
        // also do the epoch partial A*dq_dt computation in the same loop
		int state_size = 6;
		for (int i = 0; i < state_size; i++) {
            d_state[42 + i] = epoch_partial_sum_term[i];
        	for (int j = 0; j < state_size; j++) {
                unsigned int current_index = state_size + i * state_size + j;

                d_state[42 + i] += A[i * 6 + j] * state[42 + j];
	            d_state[current_index] = 0.0;
        	    for (int k = 0; k < state_size; k++) {
                    d_state[current_index] += A[i * state_size + k] * state[state_size + state_size * k + j];
                    }
                }
            }
    }
};


// Primary Acceleration
void NBODY::get_primary_acceleration(const state_type &state, state_type &acc, state_type &A_subset){
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
	
	if (prop_STM){	
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
	}
};

// Perturbing acceleration
void NBODY::get_perturbing_acceleration(const double t, SpiceBody pert_body,
	int central_body_ID, cr3bp_system sys, const state_type &state,
	state_type &acc, state_type &A_subset, state_type &sum_term){
	
    // Nomenclature: i is central body, j is perturbing body, s is spacecraft
    // Vector definition: rij = rj-ri

	// variables
	const double mass (pert_body.mu/G/sys.m_star);
	const double t_star (sys.t_star);
	SpiceDouble r_i_j [6], lighttimes;
	SpiceDouble epoch_dim (t*t_star);
	
	// get position of perturbing body WRT central body	
	spkezr_c(to_string(pert_body.ID).c_str(), epoch_dim, "J2000", "NONE",
		 to_string(central_body_ID).c_str(), r_i_j, &lighttimes);

	// non dimensionalize planet position
	r_i_j[0] = r_i_j[0]/sys.l_star;
	r_i_j[1] = r_i_j[1]/sys.l_star;
	r_i_j[2] = r_i_j[2]/sys.l_star;
    r_i_j[3] = r_i_j[3]/sys.l_star * sys.t_star;
    r_i_j[4] = r_i_j[4]/sys.l_star * sys.t_star;
    r_i_j[5] = r_i_j[5]/sys.l_star * sys.t_star;
	
	state_type r_j_s;
	double r_j_s_n, r_j_s_n_3, r_j_s_n_5, x, y, z;
	double x_o, y_o, z_o;
	double r_i_j_n, r_i_j_n_3;
	
	// Vector pointing from perturbing body to spacecraft
	r_j_s = {state[0] - r_i_j[0], state[1] - r_i_j[1], state[2] - r_i_j[2]};


	// norm of spacecraft to central body vector
	r_j_s_n = sqrt(pow(r_j_s[0], 2) + pow(r_j_s[1], 2)+ pow(r_j_s[2], 2));
	r_j_s_n_3 = pow(r_j_s_n, 3);
	r_j_s_n_5 = pow(r_j_s_n, 5);

	// norm of perturbing to central body vector
	r_i_j_n = sqrt(pow(r_i_j[0], 2) + pow(r_i_j[1], 2)+ pow(r_i_j[2], 2));
	r_i_j_n_3 = pow(r_i_j_n,3);

	// states of r_j_s vector
	x = r_j_s[0];
	y = r_j_s[1];
	z = r_j_s[2];

	// states of r_i_j vector
	x_o = r_i_j[0];
	y_o = r_i_j[1];
	z_o = r_i_j[2];

	// add perturbing accelerations to acceleration vector
	acc[0] -= mass * (x / r_j_s_n_3 + x_o / r_i_j_n_3);
	acc[1] -= mass * (y / r_j_s_n_3 + y_o / r_i_j_n_3);
	acc[2] -= mass * (z / r_j_s_n_3 + z_o / r_i_j_n_3);
	
	if (prop_STM){
	// Add perturbing contributions to A matrix
		A_subset[0] += mass * (3 * pow(x, 2) / r_j_s_n_5 - 1 / r_j_s_n_3);
		A_subset[1] += mass * 3 * x * y / r_j_s_n_5;
		A_subset[2] += mass * 3 * x * z /r_j_s_n_5;
		A_subset[3] += mass * 3 * x * y / r_j_s_n_5;
		A_subset[4] += mass * (3 * pow(y, 2) / r_j_s_n_5 - 1 / r_j_s_n_3);
		A_subset[5] += mass * (3 * y * z / r_j_s_n_5);
		A_subset[6] += mass * 3 * x * z /r_j_s_n_5;
		A_subset[7] += mass * (3 * y * z / r_j_s_n_5);
		A_subset[8] += mass * (3 * pow(z, 2) / r_j_s_n_5 - 1 / r_j_s_n_3);

    // epoch summation contribution
    state_type RIJ = {r_i_j[0], r_i_j[1], r_i_j[2], r_i_j[3], r_i_j[4], r_i_j[5]};
    states_wrt_epoch_summation(r_j_s, RIJ, mass, sum_term);
	}
};

void NBODY::states_wrt_epoch_summation(state_type RJS, state_type RIJ,
                                       const double mass, state_type &sum_term){
    // I = central body, J = Perturbing body, S = Spacecraft
    // Initialize matrix of partials del dq / delRIJ 
    std::vector<double> PM(9, 0);
    double RJS_3_inv, RJS_5_inv, RIJ_3_inv, RIJ_5_inv;

    // create variables for terms that are used frequently
    RJS_3_inv = pow(RJS[0]*RJS[0]+ RJS[1]*RJS[1]+ RJS[2]*RJS[2], -3.0/2.0);
    RJS_5_inv = pow(RJS[0]*RJS[0]+ RJS[1]*RJS[1]+ RJS[2]*RJS[2], -5.0/2.0);
    RIJ_3_inv = pow(RIJ[0]*RIJ[0]+ RIJ[1]*RIJ[1]+ RIJ[2]*RIJ[2], -3.0/2.0);
    RIJ_5_inv = pow(RIJ[0]*RIJ[0]+ RIJ[1]*RIJ[1]+ RIJ[2]*RIJ[2], -5.0/2.0);

    // fill in the partial matrix PM.
    PM[0] = -1 * mass * (-1 * RJS_3_inv + 3 * RJS[0] * RJS[0] * RJS_5_inv + RIJ_3_inv - 3 * RIJ[0] * RIJ[0] * RIJ_5_inv);
    PM[1] = -3 * mass * (RJS[0] * RJS[1] * RJS_5_inv - RIJ[0] * RIJ[1] * RIJ_5_inv);
    PM[2] = -3 * mass * (RJS[0] * RJS[2] * RJS_5_inv - RIJ[0] * RIJ[2] * RIJ_5_inv);
    PM[3] = PM[1];
    PM[4] = -1 * mass * (-1 * RJS_3_inv + 3 * RJS[1] * RJS[1] * RJS_5_inv + RIJ_3_inv - 3 * RIJ[1] * RIJ[1] * RIJ_5_inv);
    PM[5] = -3 * mass * (RJS[1] * RJS[2] * RJS_5_inv - RIJ[1] * RIJ[2] * RIJ_5_inv);
    PM[6] = PM[2];
    PM[7] = PM[5];
    PM[8] = -1 * mass * (-1 * RJS_3_inv + 3 * RJS[2] * RJS[2] * RJS_5_inv + RIJ_3_inv - 3 * RIJ[2] * RIJ[2] * RIJ_5_inv);
    std::vector<double> sum_temp(6,0);

    // partial matrix * v_ci in elementwise calculation
    sum_temp[3] = RIJ[3] * PM[0] + RIJ[4] * PM[1] + RIJ[5] * PM[2];
    sum_temp[4] = RIJ[3] * PM[3] + RIJ[4] * PM[4] + RIJ[5] * PM[5];
    sum_temp[5] = RIJ[3] * PM[6] + RIJ[4] * PM[7] + RIJ[5] * PM[8];

    sum_term[0] += sum_temp[0];
    sum_term[1] += sum_temp[1];
    sum_term[2] += sum_temp[2];
    sum_term[3] += sum_temp[3];
    sum_term[4] += sum_temp[4];
    sum_term[5] += sum_temp[5];
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

void NBODY::propagate(double t_end, double step_size, double rtol, double atol, PropObserver &o){
	namespace pl = std::placeholders;
	state_type states_and_times;
	
	// check propagation direction
	if (t_end < base_epoch_nd) {step_size = - step_size;}

	// define propagator scheme
	typedef runge_kutta_fehlberg78<state_type> rk78;

	// Create Stepper
	auto stepper = make_controlled<rk78>(atol, rtol);

	// Integrate the EOMs, capture steps with observer function o
	size_t steps;
	steps = integrate_adaptive(stepper, std::bind(&NBODY::EOM_STM, *this, pl::_1, pl::_2, pl::_3),
		       IC, base_epoch_nd, t_end, step_size, std::ref(o));
};

class PYNBODY{
	// Class that is interfaced with by Python
	public:
	// Constructor
	PYNBODY(state_type &IC, SpiceEpoch epoch, const double m, const double l, double TOF){
		double G;
		G = (6.67430e-11/pow(1000,3)); //  # km3/kg/s2	
		// Define bodies
		SpiceBody Earth("EARTH", 399, 3.9860043543609593e+05);
		SpiceBody Moon("MOON", 301, 4.9028000661637961e+03);
		SpiceBody Sun("SUN", 10, 1.3271244004193930e+11);
		SpiceBody Jupiter("JUPITER BARYCENTER", 5, 1.2671276480000021e+08);

		// Initialize our N-body system
		NBODY nbody_sys(IC, epoch, m, l, Earth, {Moon, Sun, Jupiter});
		
		// Create observer class
		PropObserver observer{}; 

		// Run the propagator
		nbody_sys.propagate(nbody_sys.base_epoch_nd + TOF, 1e-3, 1e-12, 1e-12, observer);

		// Save states and times to class attributes
		x_states = observer.x;
		t_states = observer.t;
	}

	// Attributes
	state_type t_states;
	std::vector<std::vector<double>> x_states;
};

#ifdef PYTHON_COMPILE
PYBIND11_MODULE(casper, m) {
    m.doc() = "CASPER - Python ephemeris model integrated in C++"; // optional module docstring

  //  py::class_<SpiceBody>(m, "SpiceBody")
//	    .def(py::init<const string &, const int &, const double &>())
//	    .def_readonly("name", &SpiceBody::name)
//	    .def_readonly("ID", &SpiceBody::ID)
//	    .def_readonly("mu", &SpiceBody::mu);

    //py::class_<NBODY>(m, "nbody")
//	    .def(py::init<state_type &, string, const double, const double, SpiceBody, std::vector<SpiceBody>>())
//	    .def_readonly("base_epoch", &NBODY::base_epoch_nd)
//	    .def_readonly("IC", &NBODY::IC)
//	    .def("propagate", &NBODY::propagate);

  //  py::class_<PropObserver>(m, "PropObserver")
//	    .def(py::init<>())
//	    .def_readonly("t", &PropObserver::t)
//	    .def_readonly("x", &PropObserver::x);

    py::class_<SpiceEpoch>(m, "SpiceEpoch")
        .def(py::init<SpiceDouble>())
        .def(py::init<string>());

    py::class_<PYNBODY>(m, "PyNbody")
	    .def(py::init<state_type &, SpiceEpoch, const double, const double, double>())
	    .def_readonly("x_states", &PYNBODY::x_states)
	    .def_readonly("t_states", &PYNBODY::t_states);

}
#endif

#ifndef PYTHON_COMPILE

int main(){
	cout<<"testing part now running"<<endl;

	const double G ((6.67430e-11/pow(1000,3))); //  # km3/kg/s2

	SpiceBody Earth("EARTH", 399, 3.9860043543609593e+05);
	SpiceBody Moon("MOON", 301, 4.9028000661637961e+03);
	SpiceBody Sun("SUN", 10, 1.3271244004193930e+11);
	SpiceBody Jupiter("JUPITER BARYCENTER", 5, 126712767.8578);
	// Define 48-state vector. Arbitrary Lyapunov initial state
	state_type IC_vector {0.18889952, -0.86798499, -0.34129653,  0.48008814,  0.11764799,  0.00512411,
          1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	const double m_star ((Earth.mu+Moon.mu)/G);
	const double l_star (3.8474799197904585e+05);
    SpiceEpoch epoch("March 1, 2000, 00:00:00.0000");
    cout<<epoch.utc_epoch<<endl;
	NBODY nbody_system(IC_vector, epoch, m_star, l_star, Earth, {Moon});
	cout <<"Length of IC: "<< nbody_system.n_IC << endl;
	nbody_system.EOM_STM(nbody_system.IC, nbody_system.IC, nbody_system.base_epoch_nd);	
	PropObserver obs{};
	nbody_system.propagate(nbody_system.base_epoch_nd + 5, 1e-3, 1e-12, 1e-12, obs);
};
#endif
