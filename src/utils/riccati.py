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

from scipy.special import riccati_jn, riccati_yn

def riccati_j(ell, z):
  return riccati_jn(ell, z)[0][-1]

def riccati_n(ell, z):
  return -riccati_yn(ell, z)[0][-1]
