# This code is made available as part of the FRIB-TA Summer School:
# "A Practical Walk Through Formal Scattering Theory: Connecting Bound States,
# Resonances, and Scattering States in Exotic Nuclei and Beyond," held
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

class G_0:
  def __init__(self, sys):
    self.mass = sys.mass
    self.mu = sys.mu

  def __call__(self, E, q):
    return self.mass / (self.mass * E - q**2)

  def residue(self, q0):
    return -self.mu / q0
