# As promised in class yesterday, I am providing here some numbers that you can use to check/benchmark your final-project implementations.  
# Specifically, I quote below the ground state energies E_0 for the same attractive Gaussian potential that was used in previous homework assignments, with 
# parameters V_0 = -4.0 and range R = 2.0.

# Note that even for two particles the ground-state energy below differs from the homework result because that calculation was performed for the 3D system with angular momentum 0, while now we are considering the system in 1D.

# n = 2 particles, E_0 = -3.094
# n = 3 particles, E_0 = -9.738
# n = 4 particles, E_0 = -20.046

# my results
# n=2 E_0: -3.3374750821993797
# n=3 E_0: -10.205184223506715
# n=4 E_0: -20.66960824506848
# n=5 E_0: -34.84654721518403
# (qmc) [src] python qmc_cli.py -n 6 -m 1000 -x 3000
# n=6 E_0: -52.67129155385827