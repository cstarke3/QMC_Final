# primary entry-point to run simulations.
# will handle options passed in from the command line

from simple_model import QMC
from utils.potential import V_Gauss
import argparse
import sys

# globals (yeah, I know)
particles = 2          # number of particles to simulate
dim = 1                # number of dimensions (D=1 for this exercise)
min_replicas = 500     # minimum number of replicas
max_replicas = 2000    # maximum number of replicas
max_steps = 10       # maximum number of time steps to run the simulation (τ0 = 1000)
delta_t = 0.1        # time step size (Δτ = 0.1)
xmin = -20           # minimum value of the spatial coordinate (xmin = −20)
xmax = 20            # maximum value of the spatial coordinate (xmax = 20)
bins = 200           # number of spatial bins for sorting the replicas (only used during 'Counting' to plot the ground state wave function)
seed = 42            # seed value for the random number generators (for repeatability)
mass = 1             # mass of the particle (m = 1)
E_ref = 0            # reference energy (E_ref = 0)

steps = 1000

def run_simulation():
    V_0 = -4.0
    R = 2.0
    qmc = QMC() # create a QMC object, all default values apply
    qmc.V = V_Gauss(sys, V_0, R)
    qmc.min_replicas = min_replicas
    qmc.n_part = particles
    qmc.initialize() # initialize the simulation
    for i in range(max_steps):
        V_avg, N, N_prev = qmc.step() # run the simulation forward one step
        print(f"Timestep {i} V_avg: {V_avg} N: {N}  N_prev: {N_prev}")
        
    # Energy = val.Ground_State_Energy()
    # print(Energy)

parser = argparse.ArgumentParser(
                    prog='qmc_cli.py',
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                    description="""
This script is a command line interpreter to run the QMC/DMC algorithm as part of the PY 525 fall 2023 class final.

Example usage:
    python qmc_cli.py               ::: print the help message and exit
    >>> 

""", 
                    # epilog='----'
                    )
parser.add_argument('-n', '--particles', required=True, help='the number of particles to simulate')
parser.add_argument('-m', '--min_replicas', help='the minimum number of replicas to use during the simulation')
parser.add_argument('-x', '--max_replicas', action='store_true', help='the maximum number of replicas to use during the simulation')
parser.add_argument('-s', '--steps',  help='the number of timesteps to use during the simulation')

parser.add_argument('-l', '--xmin', action='store_true', help='the left-most boundary on the x-axis')
parser.add_argument('-r', '--xmax', action='store_true', help='the right-most boundary on the x-axis')
parser.add_argument('-b', '--bins', action='store_true', help='the number of spatial “boxes” (nb) for sorting the replicas during their sampling.')


if __name__ == "__main__":

    args = parser.parse_args()

    if args.particles is not None:
        print (f"args.particles: {args.particles}")
        particles = int(args.particles)
    
    if args.steps is not None:
        print (f"args.steps: {args.steps}")
        max_steps = int(args.steps)
        
    if args.min_replicas is not None:
        print (f"args.min_replicas: {args.min_replicas}")
        min_replicas = int(args.min_replicas)
        
    run_simulation()
    

