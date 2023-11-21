# primary entry-point to run simulations.
# will handle options passed in from the command line

from model import QMC
from utils import potential
import argparse
import sys



def main():
    QMC(V, 
        particles, 
        min_replicas=500, 
        max_replicas=2000, 
        max_steps=1000, 
        delta_t=0.1, 
        xmin=-20, 
        xmax=20, 
        bins=200)

    val = QMC(potential, 2, 10)
    Energy = val.Ground_State_Energy()
    print(Energy)

    

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
parser.add_argument('-m', '--min_replicas', action='store_true', help='the minimum number of replicas to use during the simulation')
parser.add_argument('-x', '--max_replicas', action='store_true', help='the maximum number of replicas to use during the simulation')
parser.add_argument('-s', '--steps', action='store_true', help='the number of timesteps to use during the simulation')

parser.add_argument('-l', '--xmin', action='store_true', help='the left-most boundary on the x-axis')
parser.add_argument('-r', '--xmax', action='store_true', help='the right-most boundary on the x-axis')
parser.add_argument('-b', '--bins', action='store_true', help='the number of spatial “boxes” (nb) for sorting the replicas during their sampling.')




if __name__ == "__main__":

    args = parser.parse_args()
    # print(f"args.specdir: {args.specdir} args.list: {args.list}  args.pickle: {args.pickle}  args.all: {args.all}")

    if args.particles is not None:
        print (f"args.particles: {args.particles}")
    
    if args.max_replicas is not None:
        print (f"args.max_replicas: {args.max_replicas}")
        
    parser.print_help(sys.stderr)
    sys.exit(1)
    

