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

#ifdef PYTHON_COMPILE
namespace py = pybind11;
#endif

typedef std::vector<double> state_type;

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
		state_type IC; // initial condition, 42-vector
		state_type IC_nd; // initial condition in non-dim coordinates
		SpiceDouble base_epoch; // base epoch, in seconds ephemeris time
		SpiceDouble base_epoch_nd; // base epoch, in non-dim seconds
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
};

// Define NBODY class Constructor
NBODY::NBODY(state_type &ic, string str_epoch, const double m, const double l,
	SpiceBody central, std::vector<SpiceBody> perturbing):
	IC(ic),
	sys(m, l),
	central_body(central),
	perturbing_bodies(perturbing)
	{
	base_epoch = str2et(str_epoch);
	// non-dimensionalize the epoch and initial conditions
	base_epoch_nd = base_epoch/sys.t_star;
	IC_nd = IC;
        IC_nd[0] = IC_nd[0]/sys.l_star;
	IC_nd[1] = IC_nd[1]/sys.l_star;
        IC_nd[2] = IC_nd[2]/sys.l_star;
	IC_nd[3] = IC_nd[3]/sys.l_star*sys.t_star;
	IC_nd[4] = IC_nd[4]/sys.l_star*sys.t_star;
       	IC_nd[5] = IC_nd[5]/sys.l_star*sys.t_star;
	load_kernels(); // make sure kernels are loaded
	}

// EOMs with STM propagation
void NBODY::EOM_STM(state_type &state, state_type &d_state, const double t){
	state_type acc(3,0);
  	state_type A_subset (9,0);
	
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
	print_vector(d_state);
      };

// Primary Acceleration
void NBODY::get_primary_acceleration(const state_type &state, state_type &acc, state_type &A_subset){
	const double G ((6.67430e-11/pow(1000,3)));
	const double mass (central_body.mu/G/sys.m_star);
	state_type r_sc_b {state[0], state[1], state[2]};

	const double r_sc_b_n = sqrt(r_sc_b[0]*r_sc_b[0] + r_sc_b[1]*r_sc_b[1] + r_sc_b[2]*r_sc_b[2]);
	const double r_sc_b_n_5 = pow(r_sc_b_n, 5);
	const double r_sc_b_n_3 = pow(r_sc_b_n, 3);

	double x, y, z;
	x = r_sc_b[0];
	y = r_sc_b[1];
	z = r_sc_b[2];

	acc[0] -= mass * x / r_sc_b_n_3;
	acc[1] -= mass * y / r_sc_b_n_3;
	acc[2] -= mass * z / r_sc_b_n_3;

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

	const double G ((6.67430e-11/pow(1000,3)));
	const double mass (pert_body.mu/G/sys.m_star);
	const double t_star (sys.t_star);
	SpiceDouble body_position [3], lighttimes;
	SpiceDouble epoch_dim (t*t_star);
	
	spkpos_c(to_string(pert_body.ID).c_str(), epoch_dim, "J2000", "NONE", to_string(central_body_ID).c_str(), body_position, &lighttimes);
	//cout<< body_position[0]<<endl;
	body_position[0] = body_position[0]/sys.l_star;
	body_position[1] = body_position[1]/sys.l_star;
	body_position[2] = body_position[2]/sys.l_star;
	
	state_type r_sc_body, r_body_obs;
	double r_sc_b_n, r_sc_b_n_3, r_sc_b_n_5, x, y, z;
	double x_o, y_o, z_o;
	double r_b_o_n, r_b_o_n_3;

	r_sc_body = {state[0] - body_position[0], state[1] - body_position[1], state[2] - body_position[2]};
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
	// state of r_body_obs vector
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



#ifdef PYTHON_COMPILE
PYBIND11_MODULE(casper, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    py::class_<SpiceBody>(m, "SpiceBody")
	    .def(py::init<const string &, const int &, const double &>())
	    .def_readonly("name", &SpiceBody::name)
	    .def_readonly("ID", &SpiceBody::ID)
	    .def_readonly("mu", &SpiceBody::mu);

    m.def("str2et", &str2et);

    //py::class_<NBODY>(m, "nbody")
//	    .def(py::init<state_type, string >())
//	    .def_readonly("base_epoch", &NBODY::base_epoch);

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
	state_type IC_vector {1, 2, 3, 4, 5, 6, 1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1};
	const double m_star ((Earth.mu+Moon.mu)/G);
	const double l_star (3.8474799197904585e+05);
	NBODY nbody_system(IC_vector, "May 2, 2022", m_star, l_star, Earth, {Moon});
	nbody_system.EOM_STM(nbody_system.IC_nd, nbody_system.IC_nd, nbody_system.base_epoch_nd);	

};
