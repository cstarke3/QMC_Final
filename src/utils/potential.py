# This code is made available as part of the FRIB-TA Summer School:
# A Practical Walk Through Formal Scattering Theory: Connecting Bound States,
# Resonances, and Scattering States in Exotic Nuclei and Beyond, held
# virtually August 4-6, 2021.
#
# https://fribtascattering.github.io/
#
# Organizers/Lecturers:
# - Kévin Fossez (MSU/ANL)
# - Sebastian König (NCSU)
# - Heiko Hergert (MSU)
#
# Author: Sebastian König
#
# Last update: Aug 11, 2021

from .system import System
from .riccati import riccati_j

import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import quad

class Potential:
  def __init__(self, sys, name):
    self.sys = sys
    self.name = name

  def __repr__(self):
    return self.name

  def get(self, ell, p, q):
    raise RuntimeError("Abstract method call!")

  def show(self, rep='p', ell=0, **kwargs):
    if rep != 'p':
      raise RuntimeError("Unsupported potential representation!")

    m = kwargs.get("max", 10.0)
    n = kwargs.get("num", 32)

    xs = np.linspace(1e-6, m, num=n)
    px, py = np.meshgrid(xs, xs)
    v = np.vectorize(self.get, excluded=["ell"])

    plt.xlabel("p")
    plt.ylabel("p\'")
    plt.text(0.9 * m, 0.9 * m, "ell = %d" % (ell), ha="right", va="top")
    plt.contourf(px, py, v(ell, px, py), levels=n)
    plt.axis("scaled")
    plt.title(str(self))
    plt.show()

class LocalPotential(Potential):
  def __init__(self, sys, name):
    super().__init__(sys, name)

  def __call__(self, r):
    raise RuntimeError("Abstract method call!")

  def get(self, ell, p, q):
    # For convenience we define here a generic method to calculate partial-
    # wave projected momentum-space matrix elements.  Subclasses may override
    # this implementation to improve speed and/or accuracy.

    return 4.0 * np.pi / (q * p) * quad( \
      lambda r: riccati_j(ell, p * r) * self(r) * riccati_j(ell, q * r), \
      0.0, np.inf \
    )[0]

  def show(self, rep='r', ell=0, **kwargs):
    if rep == 'r':
      m = kwargs.get("max", 10.0)
      n = kwargs.get("num", 32)

      xs = np.linspace(1e-6, m, num=n)

      plt.xlabel("r")
      plt.ylabel("V(r)")
      plt.plot(xs, self(xs))
      plt.title(str(self))
      plt.show()
    elif rep == 'p':
      super().show(rep=rep, ell=ell, **kwargs)
    else:
      raise RuntimeError("Unknown potential representation!")

class V_Zero:
  def __init__(self, sys):
    super().__init__(sys, "V_Zero")

  def __call__(self, r):
    return 0.0

class V_Gauss(LocalPotential):
  def __init__(self, sys, V0, R):
    super().__init__(sys, "V_Gauss")

    self.V0 = V0
    self.R = R

    self.R2 = R**2

  def __call__(self, r):
    return self.V0 * np.exp(-r**2 / self.R2)

  def get(self, ell, p, q):
    if ell == 0:
      A = p**2 + q**2
      B = 2 * p * q
      return self.V0 * 4.0 * pow(np.pi, 1.5) * (self.R / B) \
        * np.exp(-0.25 * A * self.R2) \
        * np.sinh(0.25 * B * self.R2)
    else:
      return super().get(ell, p, q)

  def __repr__(self):
    return "V_Gauss(V0=%s, R=%s)" % (self.V0, self.R)
