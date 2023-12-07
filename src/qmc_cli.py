# primary entry-point to run simulations.
# will handle options passed in from the command line

from simple_model import QMC
from view import plot_data, mean, plot_histogram, mean_stddev, plot_energy_vs_alpha
from utils.potential import V_Gauss
import argparse
import sys
import time


# globals (yeah, I know)
particles = 2          # number of particles to simulate
dim = 1                # number of dimensions (D=1 for this exercise)
min_replicas = 500     # minimum number of replicas
max_replicas = 3000    # maximum number of replicas
max_steps = 100        # maximum number of time steps to run the simulation (τ0 = 1000)
delta_t = 0.1          # time step size (Δτ = 0.1)
xmin = -20             # minimum value of the spatial coordinate (xmin = −20)
xmax = 20              # maximum value of the spatial coordinate (xmax = 20)
bins = 100             # number of spatial bins for sorting the replicas (only used during 'Counting' to plot the ground state wave function)
seed = 42              # seed value for the random number generators (for repeatability)
mass = 1               # mass of the particle (m = 1)
E_ref = 0              # reference energy (E_ref = 0)

steps = 1000
plot = False
DEBUG = False
global_alpha = 0.13
loop = None
search_alpha = False
early_breakout = False
potential = None
samp_pct = 0.1

def run_simulation(alpha=global_alpha):
    V_0 = -4.0
    R = 2.0
    V = V_Gauss(sys, V_0, R)

    # print(f"V_Gauss: {V(1.0)}")
    qmc = QMC(V=V, 
              min_replicas=min_replicas, 
              max_replicas=max_replicas, 
              DEBUG=DEBUG,
              particle_count=particles, 
              bins=bins, 
              seed=seed, 
              alpha=alpha) 

    E_refs = []
    N_vals = []
    eyes = range(max_steps)
    epsilon = 0.0000001
    for i in eyes:
        qmc.step() # run the simulation forward one step
        Nratio = 1-qmc.N/qmc.N_target
        if DEBUG: print(f"step: {i}  E_ref: {qmc.E_ref:.6f}  N: {qmc.N}  Nratio: {Nratio:.5f}")
        E_refs.append(qmc.E_ref)
        N_vals.append(qmc.N)
        _, stddev, range_start = mean_stddev(E_refs, samp_pct)

        if early_breakout and i > 300 and abs(Nratio) < epsilon: 
            print(f" CONDITION MET @ step: {i}  E_ref: {qmc.E_ref:.6f}  N: {qmc.N}  Nratio: {Nratio:.5f}")
            break

    E_0_mean, E_0_stddev, _ = mean_stddev(E_refs, samp_pct)
    print(f"n={qmc.particle_count} E_0: {E_0_mean:.4f} +/- {E_0_stddev:.4f}  N:{qmc.N} ")
    if plot: 
        # hist, centroids = qmc.Binning()
        centroids = None
        plot_data(E_refs, N_vals, centroids, bins, title = f"QMC: n={qmc.particle_count}  N (final)={qmc.N}  alpha={alpha}",samp_pct=samp_pct)
        
    return E_0_mean, E_0_stddev, qmc.N

def find_alpha():
    alpha_lo = 0.1
    alpha_hi = 1.0
    steps = 1000
    alpha_del = (alpha_hi - alpha_lo)/steps
    alpha_xs = []
    energy_ys = []
    print(f"alpha_lo: {alpha_lo}  alpha_hi: {alpha_hi}  alpha_del: {alpha_del}  steps: {steps}")
    for i in range(steps):
        global_alpha = alpha_lo + i * alpha_del
        E_0_mean, E_0_stddev, Nval = run_simulation(alpha=global_alpha)
        alpha_xs.append(global_alpha)
        energy_ys.append(E_0_mean)
        
    plot_energy_vs_alpha(alpha_xs, energy_ys, title="Energy vs Alpha")


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
parser.add_argument('-n', '--particles',  help=f'the number of particles to simulate (default: {particles})')
parser.add_argument('-m', '--min_replicas', help=f'the minimum number of replicas to use during the simulation (default: {min_replicas})')
parser.add_argument('-s', '--steps',  help=f'the number of timesteps to use during the simulation (default: {max_steps})')

parser.add_argument('-r', '--random', help=f'set the random seed value (default: {seed})')
parser.add_argument('-t', '--trandom', action='store_true', help='set the random seed based on the current timestamp (default: varies)')

parser.add_argument('-b', '--bins', help=f'the number of spatial “boxes” (nb) for sorting the replicas during their sampling (default: {bins})')
parser.add_argument('-p', '--plot', action='store_true', help=f'plot the data (default: {plot})')
parser.add_argument('-d', '--debug', action='store_true', help=f'print out a bunch of stuff each time through the loop (default: {DEBUG})')
parser.add_argument('-a', '--alpha', help=f'modify the rate at which N/N_0 impacts potential calculation (default: {global_alpha})')

parser.add_argument('-l', '--loop', help=f'loop trhough the algorihm for n=2-10 (default: {loop})')

parser.add_argument('-g', '--gda', action='store_true', help=f'loop through simulation across a range of alphas, to see if any gets us close to E_0 = -3.10634 +/- 0.0730 (default: {search_alpha})')
parser.add_argument('-e', '--early', action='store_true', help=f'early breakout(default: {early_breakout})')

parser.add_argument('-x', '--samp_pct', help=f'sample percentage for mean_stddev(default: {samp_pct})')

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
    
    if args.plot is not None:
        # print (f"args.plot: {args.plot}")
        plot = bool(args.plot)
        
    if args.debug is not None:
        # print (f"args.debug: {args.debug}")
        DEBUG = bool(args.debug)
        
    if args.random is not None:
        # print (f"args.random: {args.random}")
        seed = int(args.random)
        
    if args.trandom:
        seed = int(time.time())
        # print (f"trandom seed: {seed}")

    if args.alpha is not None:
        # print (f"args.alpha: {args.alpha}")
        global_alpha = float(args.alpha)
       
    if args.loop is not None:
        # print (f"args.loop: {args.loop}")
        loop = int(args.loop) 
        
    if args.gda:
        search_alpha = bool(args.gda)
        # print (f"args.gda: {args.gda}")
        
    if args.early:
        early_breakout = bool(args.early)
        # print (f"args.early: {args.early}")
        
    if args.samp_pct is not None:
        # print (f"args.alpha: {args.alpha}")
        samp_pct = float(args.samp_pct)

    if loop is not None: 
        for i in range(2, loop+1):
            particles = i
            run_simulation(alpha=global_alpha)
    elif search_alpha:
        find_alpha()
        
    else:
        run_simulation(alpha=global_alpha)

