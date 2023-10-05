import casper
import bodies as bo
import crtbp

sys = crtbp.EMsys
Earth = casper.SpiceBody("EARTH", 399, 3.9860043543609593e+05);
Moon = casper.SpiceBody("MOON", 301, 4.9028000661637961e+03);
Sun = casper.SpiceBody("SUN", 10, 1.3271244004193930e+11);
Jupiter = casper.SpiceBody("JUPITER BARYCENTER", 5, 126712767.8578);

IC = [1.05903, -0.067492, -0.103524, -0.170109, 0.0960234, -0.135279, 1,
          1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]

nbody = casper.nbody(IC, "May 2, 2022", sys.m_star, sys.l_star, Earth, [Moon])
PO = casper.PropObserver()
nbody.propagate(nbody.base_epoch + 3, 1e-5, 1e-12, 1e-12, PO)
print(PO.t)
