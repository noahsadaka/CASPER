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
#include "NBodyModel.h"

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


class NBODY{
    /*
     A Class for defining an initial state in an N-body gravitational 
     point-mass model and propagating it
     * */
	public:
		
		// Attributes
		state_type IC; // initial condition, 6 or 48-vector, non-dim
		const unsigned int n_IC; // number of states, must be 6 or 48
		bool prop_STM;
		SpiceDouble base_epoch_nd; // base epoch, in non-dim seconds ephemeris time
        NBodyModel nbd_model;
		cr3bp_system sys; // CR3BP object used for non-dimensionalization
		SpiceBody central_body; // Body at center of inertial J2000 frame
		std::vector<SpiceBody> perturbing_bodies; // Perturbing bodies
        const double G = (6.67430e-11/1e9); // G converted to km3/kg/s2
		
		// Constructor Declaration
		NBODY(state_type &ic, SpiceEpoch epoch, NBodyModel nbd);

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
		void propagate(double TOF, double step_size, double rt, double at, PropObserver &o, std::vector<double> times);
};

// Define NBODY class Constructor
NBODY::NBODY(state_type &ic, SpiceEpoch epoch, NBodyModel nbd):
	// member definition, necessary for the cr3bp class
	IC(ic),
	sys(nbd.sys),
	n_IC(IC.size()),
    nbd_model(nbd),
	central_body(nbd.central_body),
	perturbing_bodies(nbd.perturbing_bodies)
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
	state_type RIS {state[0], state[1], state[2]};

	const double RIS_n = sqrt(RIS[0]*RIS[0] + RIS[1]*RIS[1] + RIS[2]*RIS[2]);
	const double RIS_n_5 = pow(RIS_n, 5);
	const double RIS_n_3 = pow(RIS_n, 3);
	
	// makes it easier to write and read
	double x, y, z;
	x = RIS[0];
	y = RIS[1];
	z = RIS[2];
	
	// acceleration terms
	acc[0] -= mass * x / RIS_n_3;
	acc[1] -= mass * y / RIS_n_3;
	acc[2] -= mass * z / RIS_n_3;
	
	if (prop_STM){	
		// STM terms
		A_subset[0] += mass * (3 * x * x / RIS_n_5 - 1 / RIS_n_3);
		A_subset[1] += mass * 3 * x * y / RIS_n_5;
		A_subset[2] += mass * 3 * x * z / RIS_n_5;
		A_subset[3] += mass * 3 * x * y / RIS_n_5;
		A_subset[4] += mass * (3 * y * y / RIS_n_5 - 1 / RIS_n_3);
		A_subset[5] += mass * 3 * y * z / RIS_n_5;
		A_subset[6] += mass * 3 * x * z / RIS_n_5;
		A_subset[7] += mass * 3 * y * z / RIS_n_5;
		A_subset[8] += mass * (3 * z * z / RIS_n_5 - 1 / RIS_n_3);
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
	SpiceDouble epoch_dim (t * t_star);
    state_type RIJ(6, 0);
	
	// get position of perturbing body WRT central body	
	spkezr_c(to_string(pert_body.ID).c_str(), epoch_dim, "J2000", "NONE",
		 to_string(central_body_ID).c_str(), r_i_j, &lighttimes);

	// non dimensionalize planet position
	RIJ[0] = r_i_j[0] / sys.l_star;
	RIJ[1] = r_i_j[1] / sys.l_star;
	RIJ[2] = r_i_j[2] / sys.l_star;
    RIJ[3] = r_i_j[3] / sys.l_star * sys.t_star;
    RIJ[4] = r_i_j[4] / sys.l_star * sys.t_star;
    RIJ[5] = r_i_j[5] / sys.l_star * sys.t_star;
	
	state_type RJS;
	double RJS_n, RJS_n_3, RJS_n_5, x, y, z;
	double x_o, y_o, z_o;
	double RIJ_n, RIJ_n_3;
	
	// Vector pointing from perturbing body to spacecraft
	RJS = {state[0] - RIJ[0], state[1] - RIJ[1], state[2] - RIJ[2]};

	// norm of spacecraft to central body vector
    // _n on name means norm, _n_X means norm to the power of X
	RJS_n = sqrt(RJS[0] * RJS[0] + RJS[1] * RJS[1] + RJS[2] * RJS[2]);
	RJS_n_3 = pow(RJS_n, 3);
	RJS_n_5 = pow(RJS_n, 5);

	// norm of perturbing to central body vector
	RIJ_n = sqrt(RIJ[0] * RIJ[0] + RIJ[1] * RIJ[1] + RIJ[2] * RIJ[2]);
	RIJ_n_3 = pow(RIJ_n, 3);

	// states of RJS vector
	x = RJS[0];
	y = RJS[1];
	z = RJS[2];

	// add perturbing accelerations to acceleration vector
	acc[0] -= mass * (x / RJS_n_3 + RIJ[0] / RIJ_n_3);
	acc[1] -= mass * (y / RJS_n_3 + RIJ[1] / RIJ_n_3);
	acc[2] -= mass * (z / RJS_n_3 + RIJ[2] / RIJ_n_3);
	
	if (prop_STM){
	// Add perturbing contributions to A matrix
		A_subset[0] += mass * (3 * x * x / RJS_n_5 - 1 / RJS_n_3);
		A_subset[1] += mass * 3 * x * y / RJS_n_5;
		A_subset[2] += mass * 3 * x * z / RJS_n_5;
		A_subset[3] += mass * 3 * x * y / RJS_n_5;
		A_subset[4] += mass * (3 * y * y / RJS_n_5 - 1 / RJS_n_3);
		A_subset[5] += mass * (3 * y * z / RJS_n_5);
		A_subset[6] += mass * 3 * x * z / RJS_n_5;
		A_subset[7] += mass * (3 * y * z / RJS_n_5);
		A_subset[8] += mass * (3 * z * z / RJS_n_5 - 1 / RJS_n_3);

    // epoch summation contribution
    states_wrt_epoch_summation(RJS, RIJ, mass, sum_term);
	}
};

void NBODY::states_wrt_epoch_summation(state_type RJS, state_type RIJ,
                                       const double mass, state_type &sum_term){
    // I = central body, J = Perturbing body, S = Spacecraft
    // Initialize matrix of partials del dq / delRIJ 
    std::vector<double> PM(9, 0);
    double RJS_3_inv, RJS_5_inv, RIJ_3_inv, RIJ_5_inv;

    // create variables for terms that are used frequently
    RJS_3_inv = pow(RJS[0] * RJS[0] + RJS[1] * RJS[1] + RJS[2] * RJS[2], -3.0 / 2.0);
    RJS_5_inv = pow(RJS[0] * RJS[0] + RJS[1] * RJS[1] + RJS[2] * RJS[2], -5.0 / 2.0);
    RIJ_3_inv = pow(RIJ[0] * RIJ[0] + RIJ[1] * RIJ[1] + RIJ[2] * RIJ[2], -3.0 / 2.0);
    RIJ_5_inv = pow(RIJ[0] * RIJ[0] + RIJ[1] * RIJ[1] + RIJ[2] * RIJ[2], -5.0 / 2.0);

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
    std::vector<double> sum_temp(6, 0);

    // partial matrix * v_ci as elementwise calculation
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
	std::vector<double> A(36, 0);
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

void NBODY::propagate(double t_end, double step_size, double rtol, double atol, PropObserver &o, std::vector<double> times){
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

    if ( times.size() > 0 ) {
     	steps = integrate_times(stepper, std::bind(&NBODY::EOM_STM, *this, pl::_1, pl::_2, pl::_3),
		       IC, times.begin(), times.end(), step_size, std::ref(o));
    }
    else{

	steps = integrate_adaptive(stepper, std::bind(&NBODY::EOM_STM, *this, pl::_1, pl::_2, pl::_3),
		       IC, base_epoch_nd, t_end, step_size, std::ref(o));
    }
};

class PYNBODY{
	// Class that is interfaced with by Python
	public:
	// Constructor
	PYNBODY(state_type &IC, SpiceEpoch epoch, NBodyModel nbd, double TOF, std::vector<double> t_eval){
        // Initialize the state propagation class
		NBODY nbody_sys(IC, epoch, nbd);
		
		// Create observer class
		PropObserver observer{}; 

		// Run the propagator
		nbody_sys.propagate(nbody_sys.base_epoch_nd + TOF, 1e-3, 1e-12, 1e-12, observer, t_eval);

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

    py::class_<SpiceBody>(m, "SpiceBody")
	    .def(py::init<const string &, const int &, const double &>())
	    .def_readonly("name", &SpiceBody::name)
	    .def_readonly("ID", &SpiceBody::ID)
	    .def_readonly("mu", &SpiceBody::mu);

    py::class_<SpiceEpoch>(m, "SpiceEpoch")
        .def(py::init<SpiceDouble>())
        .def(py::init<string>());

    py::class_<PYNBODY>(m, "PyNbody")
	    .def(py::init<state_type &, SpiceEpoch, NBodyModel, double, std::vector<double>>())
	    .def_readonly("x_states", &PYNBODY::x_states)
	    .def_readonly("t_states", &PYNBODY::t_states);
    
    py::class_<NBodyModel>(m, "NBodyModel")
        .def(py::init<const double, const double, SpiceBody, std::vector<SpiceBody>>());

}
#endif

#ifndef PYTHON_COMPILE

int main(){
	cout<<"testing part now running"<<endl;

	const double G ((6.67430e-11/1e9)); //  # km3/kg/s2
    // Define the NBodyModel class
	SpiceBody Earth("EARTH", 399, 3.9860043543609593e+05);
	SpiceBody Moon("MOON", 301, 4.9028000661637961e+03);
	SpiceBody Sun("SUN", 10, 1.3271244004193930e+11);
	SpiceBody Jupiter("JUPITER BARYCENTER", 5, 126712767.8578);
	const double m_star ((Earth.mu+Moon.mu)/G);
	const double l_star (3.8474799197904585e+05);
    NBodyModel nbd(m_star, l_star, Earth, {Moon, Sun, Jupiter});

	// Define 48-state vector. Arbitrary Lyapunov initial state
	state_type IC_vector {0.18889952, -0.86798499, -0.34129653,  0.48008814,  0.11764799,  0.00512411,
          1.0, 0.0, 0.0, 0.0, 0.0, 0.0, // initialize STM with
          0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // and initialize epoch partial
    std::vector<double> times;
    SpiceEpoch epoch("March 1, 2000, 00:00:00.0000");
    cout<<epoch.utc_epoch<<endl;
	NBODY nbody_system(IC_vector, epoch, nbd);
	cout <<"Length of IC: "<< nbody_system.n_IC << endl;
	nbody_system.EOM_STM(nbody_system.IC, nbody_system.IC, nbody_system.base_epoch_nd);	
	PropObserver obs{};
	nbody_system.propagate(nbody_system.base_epoch_nd + 5, 1e-3, 1e-12, 1e-12, obs, times);
};
#endif
