#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <vector>
#include <string>
#include "SpiceUsr.h"
#include "pybind11/pybind11.h"
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/io.hpp>
using namespace std; 

int kernelify(){
string fullk_path;
string gen_ken= "generic_kernels.tm";
string full_ken;
fullk_path = getenv("KERNEL_PATH");
full_ken = fullk_path + gen_ken;
furnsh_c(full_ken.c_str());
int tot;
ktotal_c("ALL", &tot);
return tot;
}

int add_numbers(int a, int b)
{
	int c;
	c = kernelify();
	return c;
}



//int main()
//{
//SpiceDouble et;
//string fullk_path;
//string gen_ken= "generic_kernels.tm";
//string full_ken;
//fullk_path = getenv("KERNEL_PATH");

//full_ken = fullk_path + gen_ken;
//furnsh_c(full_ken.c_str());
//str2et_c("September 17, 1996", &et);
//SpiceInt kernelcount;
//ktotal_c("ALL", &kernelcount);

//cout << et << endl;
//cout << kernelcount << " kernels loaded" << endl;

//int a, b, c;
//a = 1;
//b = 2;
//c = add_numbers(a,b);
//cout << c << endl;
//} 

PYBIND11_MODULE(spicetest, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("add_numbers", &add_numbers, "A function that adds two numbers");
}

