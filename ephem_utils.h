#ifndef CASPER_EPHEM_UTILS_H
#define CASPER_EPHEM_UTILS_H

// Standard Imports
#include <iostream>
#include <string>
#include <cstring>

// SPICE toolkit
#include "SpiceUsr.h"


int load_kernels(){
	/** 
	 Load SPICE kernels. Need an environment variable called KERNEL_PATH 
	 for the folder where the kernels are stored

	 Checks to see if kernels have been loaded. If not, loads them. If so,
	 does not load them

	 Returns the number of loaded kernels
	 **/
	int n_kernels;
	ktotal_c("ALL", &n_kernels);
	if (n_kernels == 0){
		string full_kernel_path, kernel_folder_path;
		string kernel_name = "generic_kernels.tm";
		kernel_folder_path = getenv("KERNEL_PATH");
		full_kernel_path = kernel_folder_path + kernel_name;
		furnsh_c(full_kernel_path.c_str());
		ktotal_c("ALL", &n_kernels);
	}
	return n_kernels;
}

double str2et(string &epoch_string){
	/**
 	Convert string to ephemeris time, accessible from Python
	**/
	load_kernels();

	SpiceDouble et;
	str2et_c(epoch_string.c_str(), &et);
	return et;
};

#endif
