import casper
import numpy as np
import bodies as bo
from crtbp import EMsys
import matplotlib.pyplot as plt
import time
import ephemeristools as et

sys = EMsys

IC = [ 0.18889952, -0.86798499, -0.34129653,  0.48008814,  0.11764799,  0.00512411,
          1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
start_time = time.time()
nbody = casper.PyNbody(IC, "March 1, 2000, 00:00:00.0000", sys.m_star, sys.l_star, 5)
print(time.time()-start_time)



t = nbody.t_states;
print(t[-1])
x = np.asarray(nbody.x_states);
plt.plot(x[:,0], x[:,1])
plt.gca().set_aspect('equal')
plt.show()
#print(nbody.t_states)
#print(nbody.x_states)
