# CASPER
C++ Aided, SPICE-based, Propagation for Ephemeris model Results (name WIP).

A C++ N-body propagation tool that is wrapped with pybind11 to be called from Python.

Uses Boost.Odeint RK78 integrator 

Very much a work in progress

## TODO
 - [x] choose whether to propagate the STM or not depending on the size of the initial condition vector
 - [x] throw exceptions if the initial condition vector is not a 6-vector or a 42-vector
 - [ ] implement event finding
 - [ ] add options to specify what the central body and perturbing bodies are from Python
 - [ ] Use fewer time steps if integrating for corrections vs integrating for detailed state history and/or plotting
 - [x] Implement some sort of a smart Spice time class so I can pass it an epoch string or an ET

## Notes to make this work
 - Need to have a folder with the meta-kernel saved in it with the name "generic_kernels.tm" assigned to an environment variable "KERNEL_PATH"
 - Need to have your python installation path added to the "CPLUS_INCLUDE_PATH" environment variable
 - boost needs to be installed
