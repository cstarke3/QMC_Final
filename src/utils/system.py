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

import numpy as np

class System:
  def __init__(self, mass=1.0, scale=1.0):
    self.mass = mass
    self.scale = scale

    self.mu = self.mass / 2.0

  def e_from_k(self, k):
    return 0.5 * k**2 / self.mu

  def k_from_e(self, e):
    return np.sqrt(2.0 * self.mu * e)

  def __repr__(self):
    return "System(mass=%s, scale=%s)" % (self.mass, self.scale)
