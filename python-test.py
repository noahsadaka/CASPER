import casper
import spiceypy as sp
from scipy.integrate import solve_ivp
import numpy as np
import bodies as bo
from crtbp import EMsys
import matplotlib.pyplot as plt
import time
import ephemeristools as et

sys = EMsys

IC = [ 0.18889952, -0.86798499, -0.34129653,  0.48008814,  0.11764799,  0.00512411,
          1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
start_time = time.time()
base_epoch = 'March 1, 2000, 00:00:00.0000'
nbody = casper.PyNbody(IC, "March 1, 2000, 00:00:00.0000", sys.m_star, sys.l_star, 5)
print(time.time()-start_time)
et.furnsh_generic_kernels()
base_epoch_ET = sp.str2et(base_epoch)
start_time = time.time()
traj = solve_ivp(et.Nbody_prop,
                 (base_epoch_ET/sys.t_star, (base_epoch_ET)/sys.t_star+5),
                 IC,
                  atol=1e-12,
                  rtol=1e-12,
                  method='RK45',
                  args=(bo.earth, [bo.moon, bo.sun, bo.jupiter_bary], sys))
print(time.time()-start_time)
t = nbody.t_states;
#print(t[-1])
x = np.asarray(nbody.x_states);
plt.figure(1)
plt.plot(x[:,0], x[:,1], label='C++')
plt.plot(traj.y[0], traj.y[1], label='Py')
plt.legend()
plt.gca().set_aspect('equal')

plt.figure(2)
plt.plot(t, x[:,0], label='C++')
plt.plot(traj.t, traj.y[0], label='Python')
plt.legend()
plt.show()




