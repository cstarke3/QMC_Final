# primary entry-point to run simulations.
# will handle options passed in from the command line

from simple_model import QMC
from view import plot_data, mean
from utils.potential import V_Gauss
import argparse
import sys
import time


# globals (yeah, I know)
particles = 2          # number of particles to simulate
dim = 1                # number of dimensions (D=1 for this exercise)
min_replicas = 500     # minimum number of replicas
max_replicas = 2000    # maximum number of replicas
max_steps = 1000       # maximum number of time steps to run the simulation (τ0 = 1000)
delta_t = 0.1          # time step size (Δτ = 0.1)
xmin = -20             # minimum value of the spatial coordinate (xmin = −20)
xmax = 20              # maximum value of the spatial coordinate (xmax = 20)
bins = 200             # number of spatial bins for sorting the replicas (only used during 'Counting' to plot the ground state wave function)
seed = 42              # seed value for the random number generators (for repeatability)
mass = 1               # mass of the particle (m = 1)
E_ref = 0              # reference energy (E_ref = 0)

steps = 1000
plot_it = False
DEBUG = False

def run_simulation():
    V_0 = -4.0
    R = 2.0
    V = V_Gauss(sys, V_0, R)
    qmc = QMC(V=V,min_replicas=min_replicas, n_part=particles, bins=bins, seed=seed) 

    E_refs = []
    eyes = range(max_steps)
    for i in eyes:
        qmc.step() # run the simulation forward one step
        if DEBUG: print(f"step: {i} E_ref: {qmc.E_ref} N: {qmc.N}")
        E_refs.append(qmc.E_ref)

    E_0 = mean(E_refs)
    print(f"n={qmc.n_part} E_0: {E_0}")
    if plot_it: plot_data(eyes,E_refs, title = f"QMC Simulation for n={qmc.n_part} steps={max_steps}")
    

parser = argparse.ArgumentParser(
                    prog='qmc_cli.py',
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                    description="""
This script is a command line interpreter to run the QMC/DMC algorithm as part of the Fall2023 PY 525 final project.

Example usage:
    python qmc_cli.py -h            ::: print the help message and exit
    python qmc_cli.py -n 6

""", 
                    # epilog='----'
                    )
parser.add_argument('-n', '--particles', required=True, help=f'the number of particles to simulate (default: {particles})')
parser.add_argument('-m', '--min_replicas', help=f'the minimum number of replicas to use during the simulation (default: {min_replicas}))')
parser.add_argument('-x', '--max_replicas', help=f'the maximum number of replicas to use during the simulation (default: {max_replicas}))')
parser.add_argument('-s', '--steps',  help=f'the number of timesteps to use during the simulation (default: {max_steps})))')

parser.add_argument('-r', '--random', help=f'set the random seed value (default: {seed}))))')
parser.add_argument('-t', '--trandom', action='store_true', help='set the random seed based on the current timestamp (default: varies))))')

parser.add_argument('-b', '--bins', help=f'the number of spatial “boxes” (nb) for sorting the replicas during their sampling (default: {bins}))')
parser.add_argument('-p', '--plot_it', action='store_true', help=f'plot the data (default: {plot_it}))')
parser.add_argument('-d', '--debug', action='store_true', help=f'plot the data (default: {DEBUG}))')


if __name__ == "__main__":

    args = parser.parse_args()

    if args.particles is not None:
        # print (f"args.particles: {args.particles}")
        particles = int(args.particles)
    
    if args.steps is not None:
        # print (f"args.steps: {args.steps}")
        max_steps = int(args.steps)
        
    if args.min_replicas is not None:
        # print (f"args.min_replicas: {args.min_replicas}")
        min_replicas = int(args.min_replicas)
    
    if args.plot_it is not None:
        # print (f"args.plot_it: {args.plot_it}")
        plot_it = bool(args.plot_it)
        
    if args.debug is not None:
        # print (f"args.debug: {args.debug}")
        DEBUG = bool(args.debug)
        
    if args.random is not None:
        # print (f"args.random: {args.random}")
        seed = int(args.random)
    if args.trandom:
        seed = int(time.time())
        # print (f"trandom seed: {seed}")
        
    run_simulation()
    

